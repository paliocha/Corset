// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.
//
// Last modified 22 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#include <Cluster.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <ranges>
#include <omp.h>

#ifdef __SSE2__
#include <emmintrin.h>   // SSE2 — guaranteed on x86_64
#endif

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::min;
using std::ofstream;
using std::ostringstream;
using std::ios_base;


// ── Distance between two transcript clusters ────────────────────────
// Returns proportion of shared reads (1 = identical, 0 = disjoint).
// If a log-likelihood ratio test detects differential expression between
// experimental groups, the distance is forced to 0 (maximum separation).

float Cluster::get_dist(int i, int j, float lrt_softness) {
    int shared_reads   = 0;
    float total_reads_i = 0;
    float total_reads_j = 0;

    const int nsamples = Transcript::samples;
    const int ngroups  = Transcript::groups;

    // Per-sample shared read counts
    vector<int> sample_shared_read(nsamples, 0);

    for (int s = 0; s < nsamples; s++) {
        // ── Early-out: skip samples where either group has zero reads ──
        // Avoids the sorted-list intersection entirely (~30-50% of samples
        // are zero for any given transcript pair in typical RNA-Seq data).
        if (read_group_sizes[i][s] == 0 || read_group_sizes[j][s] == 0) {
            total_reads_i += read_group_sizes[i][s];
            total_reads_j += read_group_sizes[j][s];
            continue;
        }

        const vector<int> &rg_i = read_groups[s][i];
        const vector<int> &rg_j = read_groups[s][j];
        const int ni = static_cast<int>(rg_i.size());
        const int nj = static_cast<int>(rg_j.size());
        const int *a = rg_i.data();
        const int *b = rg_j.data();

        int pi = 0, pj = 0;

#ifdef __SSE2__
        // ── SIMD block intersection ────────────────────────────────────
        // Process list A in blocks of 4: for each block, compare all
        // qualifying B elements via _mm_cmpeq_epi32 (4 comparisons in
        // one instruction).  Falls through to scalar for the remainder.
        while (pi + 3 < ni && pj < nj) {
            __m128i va  = _mm_loadu_si128(reinterpret_cast<const __m128i *>(a + pi));
            int a_max = a[pi + 3];
            int a_min = a[pi];

            // Process all b[pj] that fall within [a_min, a_max]
            while (pj < nj && b[pj] <= a_max) {
                if (b[pj] < a_min) { ++pj; continue; }
                __m128i vb  = _mm_set1_epi32(b[pj]);
                __m128i eq  = _mm_cmpeq_epi32(va, vb);
                int mask = _mm_movemask_ps(_mm_castsi128_ps(eq));
                if (mask)
                    sample_shared_read[s] += get_read(b[pj])->get_weight();
                ++pj;
            }
            pi += 4;
        }
#endif

        // ── Scalar remainder ───────────────────────────────────────────
        while (pi < ni && pj < nj) {
            if      (a[pi] < b[pj]) ++pi;
            else if (a[pi] > b[pj]) ++pj;
            else {
                sample_shared_read[s] += get_read(a[pi])->get_weight();
                ++pi; ++pj;
            }
        }

        shared_reads  += sample_shared_read[s];
        total_reads_i += read_group_sizes[i][s];
        total_reads_j += read_group_sizes[j][s];
    }

    float dist = shared_reads / min(total_reads_i, total_reads_j);
    if (dist == 0) return 0;

    // Log-likelihood ratio test for differential expression
    vector<float> x(ngroups, 1.0f), y(ngroups, 1.0f);

    for (int s = 0; s < nsamples; s++) {
        int g      = sample_groups[s];
        int shared = static_cast<int>(0.5f * sample_shared_read[s]);
        x[g] += read_group_sizes[i][s] - shared;
        y[g] += read_group_sizes[j][s] - shared;
    }

    float non_null = 0, x_all = 0, y_all = 0;
    for (int g = 0; g < ngroups; g++) {
        float r = y[g] / x[g];
        non_null += y[g] * log(r * x[g]) - r * x[g];
        non_null += x[g] * log(x[g])     - x[g];
        x_all += x[g];
        y_all += y[g];
    }

    float null = 0;
    float r_all = y_all / x_all;
    for (int g = 0; g < ngroups; g++) {
        float mean_x = (x[g] + y[g]) / (1 + r_all);
        null += y[g] * log(r_all * mean_x) - r_all * mean_x;
        null += x[g] * log(mean_x)         - mean_x;
    }

    float D = 2 * non_null - 2 * null;

    // Hard cutoff (hierarchical, or Leiden without --lrt-softness)
    if (lrt_softness <= 0.0f)
        return (D > D_cut) ? 0 : dist;

    // Soft sigmoid decay for Leiden: smoothly down-weight borderline
    // differential pairs instead of a binary gate at D_cut
    if (D <= 0) return dist;
    float scale = 1.0f / (1.0f + expf(lrt_softness * (D - D_cut)));
    return dist * scale;
}


// ── Merge cluster j into cluster i ──────────────────────────────────

void Cluster::merge(int i, int j) {
    // Mark j as dead
    alive_[j] = false;

    // Merge group membership
    groups[i].insert(groups[i].end(), groups[j].begin(), groups[j].end());
    groups[j].clear();

    // Merge read-group sorted lists per sample
    for (int s = 0; s < Transcript::samples; s++) {
        vector<int> &ri = read_groups[s][i];
        vector<int> &rj = read_groups[s][j];

        vector<int> merged;
        merged.reserve(ri.size() + rj.size());

        auto it1 = ri.begin(), end1 = ri.end();
        auto it2 = rj.begin(), end2 = rj.end();
        int dup_weight = 0;

        while (it1 != end1 && it2 != end2) {
            if      (*it1 < *it2) { merged.push_back(*it1); ++it1; }
            else if (*it1 > *it2) { merged.push_back(*it2); ++it2; }
            else {
                dup_weight += get_read(*it1)->get_weight();
                merged.push_back(*it1);
                ++it1; ++it2;
            }
        }
        merged.insert(merged.end(), it1, end1);
        merged.insert(merged.end(), it2, end2);

        // Incremental size update: avoid re-scanning all reads
        read_group_sizes[i][s] += read_group_sizes[j][s] - dup_weight;

        ri = std::move(merged);
        rj.clear();
        rj.shrink_to_fit();
    }
    read_group_sizes[j].clear();

    // ── Collect live neighbors of i ∪ j (Tier 1: O(degree) not O(n)) ──
    vector<int> candidates;
    {
        // Sorted merge of adj_[i] and adj_[j], skipping dead entries and i,j
        auto &ai = adj_[i];
        auto &aj = adj_[j];
        candidates.reserve(ai.size() + aj.size());

        auto p = ai.begin(), pe = ai.end();
        auto q = aj.begin(), qe = aj.end();
        while (p != pe && q != qe) {
            if (*p < *q)      { if (alive_[*p] && *p != j) candidates.push_back(*p); ++p; }
            else if (*p > *q) { if (alive_[*q] && *q != i) candidates.push_back(*q); ++q; }
            else              { if (alive_[*p] && *p != i && *p != j) candidates.push_back(*p); ++p; ++q; }
        }
        while (p != pe) { if (alive_[*p] && *p != j) candidates.push_back(*p); ++p; }
        while (q != qe) { if (alive_[*q] && *q != i) candidates.push_back(*q); ++q; }
    }

    // ── Tier 2: Parallel distance recomputation ────────────────────────
    const int ncand = static_cast<int>(candidates.size());
    vector<unsigned char> new_dists(ncand);

    // Use nested parallelism only for big clusters (avoids overhead for tiny ones)
    #pragma omp parallel for schedule(static) if(ncand > 128)
    for (int idx = 0; idx < ncand; idx++) {
        new_dists[idx] = static_cast<unsigned char>(get_dist(i, candidates[idx]) * UCHAR_MAX);
    }

    // ── Sequential: apply distance updates, rebuild adj_[i] ────────────
    // Remove j's old entries from the distance matrix and neighbors' adj lists
    for (int k : adj_[j]) {
        if (!alive_[k]) continue;
        if (k < j) dist.remove(j, k);
        else       dist.remove(k, j);
        auto &ak = adj_[k];
        auto it = std::lower_bound(ak.begin(), ak.end(), j);
        if (it != ak.end() && *it == j) ak.erase(it);
    }
    adj_[j].clear();

    // Remove i's old entries (will be re-added below if still linked)
    for (int k : adj_[i]) {
        if (!alive_[k]) continue;
        if (k < i) dist.remove(i, k);
        else       dist.remove(k, i);
        auto &ak = adj_[k];
        auto it = std::lower_bound(ak.begin(), ak.end(), i);
        if (it != ak.end() && *it == i) ak.erase(it);
    }

    // Rebuild adj_[i] with fresh distances
    adj_[i].clear();
    adj_[i].reserve(ncand);
    for (int idx = 0; idx < ncand; idx++) {
        int k = candidates[idx];
        unsigned char d = new_dists[idx];
        if (d > 0) {
            // Store in canonical order (higher index, lower index)
            if (k > i) dist.set(k, i, d);
            else       dist.set(i, k, d);
            pq_.push({d, std::max(i, k), std::min(i, k)});
            adj_[i].push_back(k);
            // Add i to k's neighbor list (maintain sorted invariant)
            auto &ak = adj_[k];
            auto pos = std::lower_bound(ak.begin(), ak.end(), i);
            if (pos == ak.end() || *pos != i)
                ak.insert(pos, i);
        }
    }
    std::ranges::sort(adj_[i]);
}


// ── Find the closest pair via lazy-deletion max-heap ────────────────

unsigned char Cluster::find_next_pair(int &max_i, int &max_j) {
    while (!pq_.empty()) {
        HeapEntry top = pq_.top();
        pq_.pop();
        if (dist.no_link(top.i, top.j)) continue;
        if (dist.get(top.i, top.j) != top.value) continue;
        max_i = top.i;
        max_j = top.j;
        return top.value;
    }
    return 0;
}


// ── Count reads per cluster for a given sample ──────────────────────

vector<int> Cluster::get_counts(int s) {
    const int nreads   = n_reads();
    const int nclusters = static_cast<int>(read_groups[s].size());

    // Thread-safe RNG seed derived from cluster id + sample
    unsigned int seed = static_cast<unsigned int>(id_ * 17 + s * 31 + 42);

    // Build inverted index: read → list of clusters it appears in
    vector<vector<int>> read_to_clusters(nreads);
    for (int g = 0; g < nclusters; g++) {
        for (int idx : read_groups[s][g])
            read_to_clusters[idx].push_back(g);
    }

    // Deduplicate each read's cluster list
    for (int rd = 0; rd < nreads; rd++) {
        auto &v = read_to_clusters[rd];
        if (v.size() > 1) {
            std::ranges::sort(v);
            auto [first, last] = std::ranges::unique(v);
            v.erase(first, last);
        }
    }

    // Randomly assign each read to one of its clusters
    vector<int> counts(nclusters, 0);
    for (int rd = 0; rd < nreads; rd++) {
        int n_align = static_cast<int>(read_to_clusters[rd].size());
        if (n_align == 0) continue;
        int weight = get_read(rd)->get_weight();
        for (; weight > 0; weight--)
            counts[read_to_clusters[rd][rand_r(&seed) % n_align]]++;
    }
    return counts;
}


// ── Main clustering loop ────────────────────────────────────────────

void Cluster::cluster(map<float, string> &thresholds,
                      const string &method_tag) {
    if (n_trans() > 1000) {
        #pragma omp critical(print)
        cout << "cluster with " << n_trans() << " transcripts.. this might take a while" << endl;
    }

    initialise_matrix();

    auto itr_d = thresholds.begin();
    for (int n = n_trans(); n > 1; n--) {
        int i = 0, j = 0;
        float distance = 1.0f - static_cast<float>(find_next_pair(i, j)) / UCHAR_MAX;

        if (n > 1000 && n % 200 == 0) {
            #pragma omp critical(print)
            cout << "down to " << n << " clusters. dist=" << distance << endl;
        }

        while (itr_d != thresholds.end() && distance > itr_d->first) {
            output_clusters(itr_d->second, method_tag);
            itr_d++;
        }
        if (itr_d == thresholds.end() || distance == 1.0f) break;

        merge(i, j);
    }

    // Report final clustering for any remaining thresholds
    while (itr_d != thresholds.end()) {
        output_clusters(itr_d->second, method_tag);
        itr_d++;
    }
}


// ── Output cluster assignments and counts ───────────────────────────

void Cluster::output_clusters(const string &threshold,
                              const string &method_tag) {
    // Compute counts for all samples
    vector<vector<int>> counts;
    counts.reserve(Transcript::samples);
    for (int s = 0; s < Transcript::samples; s++)
        counts.push_back(get_counts(s));

    // Buffer counts output
    ostringstream counts_buf;
    int id = 0;
    int nclusters = static_cast<int>(counts[0].size());
    for (int g = 0; g < nclusters; g++) {
        if (groups[g].empty()) continue;
        counts_buf << cluster_id_prefix << get_id() << cluster_id_joiner << id;
        for (int s = 0; s < Transcript::samples; s++)
            counts_buf << "\t" << counts[s][g];
        counts_buf << "\n";
        id++;
    }

    // Buffer cluster assignments
    ostringstream cluster_buf;
    id = 0;
    for (int g = 0; g < static_cast<int>(groups.size()); g++) {
        for (int t = 0; t < static_cast<int>(groups[g].size()); t++) {
            cluster_buf << get_tran(groups[g][t])->get_name()
                        << "\t" << cluster_id_prefix << get_id()
                        << cluster_id_joiner << id << "\n";
        }
        if (!groups[g].empty()) id++;
    }

    // Write atomically under critical section (RAII handles close on scope exit)
    #pragma omp critical(file_io)
    {
        string counts_fn  = file_prefix + method_tag + string(file_counts) + threshold + string(file_ext);
        string cluster_fn = file_prefix + method_tag + string(file_clusters) + threshold + string(file_ext);

        ofstream countsFile(counts_fn, ios_base::app);
        countsFile << counts_buf.str();

        ofstream clusterFile(cluster_fn, ios_base::app);
        clusterFile << cluster_buf.str();
    }
}


// ── Set up read-group data structures ───────────────────────────────
// Shared setup used by both hierarchical clustering (initialise_matrix)
// and Leiden clustering (cluster_leiden).  Assigns positional indices,
// builds per-sample read-group lists and sizes.

void Cluster::setup_read_groups() {
    const int ntrans   = n_trans();
    const int nsamples = Transcript::samples;

    // Assign positional indices to transcripts
    for (int i = 0; i < ntrans; i++)
        get_tran(i)->pos(i);

    // Clear before resize — idempotent for "both" mode where Leiden
    // calls setup_read_groups() first and hierarchical calls it again.
    read_groups.clear();
    read_groups.resize(nsamples);
    for (int s = 0; s < nsamples; s++)
        read_groups[s].resize(ntrans);
    read_group_sizes.clear();
    read_group_sizes.resize(ntrans);
    for (int t = 0; t < ntrans; t++)
        read_group_sizes[t].resize(nsamples, 0);

    // Build read-group lists
    for (int r = 0; r < n_reads(); r++) {
        Read *read   = get_read(r);
        int   sample = read->get_sample();

        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            read_groups[sample][i].push_back(r);
            read_group_sizes[i][sample] += read->get_weight();
        }
    }

    // Each transcript gets its own initial cluster group
    groups.clear();
    groups.resize(ntrans);
    for (int n = 0; n < ntrans; n++)
        groups[n].push_back(n);
}


// ── Initialise the distance matrix ──────────────────────────────────

void Cluster::initialise_matrix() {
    const int ntrans = n_trans();

    setup_read_groups();
    dist.set_size(ntrans);

    // Collect non-zero distance pairs from shared reads
    vector<std::pair<int, int>> nz_pairs;
    for (int r = 0; r < n_reads(); r++) {
        Read *read = get_read(r);
        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            for (auto t2 = read->align_begin(); t2 != t1; ++t2) {
                int j = (*t2)->pos();
                nz_pairs.emplace_back(std::max(i, j), std::min(i, j));
            }
        }
    }

    // Deduplicate pairs to avoid redundant distance calculations
    std::ranges::sort(nz_pairs);
    auto [first, last] = std::ranges::unique(nz_pairs);
    nz_pairs.erase(first, last);

    // Initialise adjacency lists and alive flags
    adj_.resize(ntrans);
    alive_.assign(ntrans, true);

    // Compute distances and seed the max-heap + adjacency lists
    dist.reserve(nz_pairs.size());
    for (auto &[hi, lo] : nz_pairs) {
        unsigned char d = static_cast<unsigned char>(get_dist(hi, lo) * UCHAR_MAX);
        dist.set(hi, lo, d);
        if (d > 0) {
            pq_.push({d, hi, lo});
            adj_[hi].push_back(lo);
            adj_[lo].push_back(hi);
        }
    }

    // Sort adjacency lists for efficient set-union in merge()
    for (int i = 0; i < ntrans; i++)
        std::ranges::sort(adj_[i]);
}


// ── Debug: print alignment details for large clusters ───────────────

void Cluster::print_alignments() const {
    if (n_trans() < 10000) return;

    for (int r = 0; r < n_reads(); r++) {
        Read *read = get_read(r);
        cout << "cluster=" << get_id() << "\tsample=" << read->get_sample();

        vector<int> positions;
        for (auto t = read->align_begin(); t != read->align_end(); ++t)
            positions.push_back((*t)->pos());
        std::ranges::sort(positions);

        cout << "\tthere are " << positions.size() << " reads:";
        for (int pos : positions)
            cout << "\t" << pos;
        cout << endl;
    }
}
