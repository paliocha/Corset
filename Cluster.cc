// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <Cluster.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ranges>
#include <omp.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::min;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::ios_base;


// ── Distance between two transcript clusters ────────────────────────
// Returns proportion of shared reads (1 = identical, 0 = disjoint).
// If a log-likelihood ratio test detects differential expression between
// experimental groups, the distance is forced to 0 (maximum separation).

float Cluster::get_dist(int i, int j) {
    int shared_reads   = 0;
    float total_reads_i = 0;
    float total_reads_j = 0;

    const int nsamples = Transcript::samples;
    const int ngroups  = Transcript::groups;

    // Per-sample shared read counts
    vector<int> sample_shared_read(nsamples, 0);

    for (int s = 0; s < nsamples; s++) {
        vector<int> &rg_i = read_groups[s][i];
        vector<int> &rg_j = read_groups[s][j];

        auto it1 = rg_i.begin(), end1 = rg_i.end();
        auto it2 = rg_j.begin(), end2 = rg_j.end();

        while (it1 != end1 && it2 != end2) {
            if      (*it1 < *it2) ++it1;
            else if (*it1 > *it2) ++it2;
            else {
                sample_shared_read[s] += get_read(*it1)->get_weight();
                ++it1; ++it2;
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
    return (D > D_cut) ? 0 : dist;
}


// ── Merge cluster j into cluster i ──────────────────────────────────

void Cluster::merge(int i, int j) {
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

    // Update distance matrix — i absorbs j
    int ntrans = n_trans();

    for (int k = 0; k < i; k++) {
        if (dist.no_link(i, k) && ((k < j && dist.no_link(j, k)) || (j < k && dist.no_link(k, j))))
            continue;
        if (dist.get(i, k) == UCHAR_MAX &&
            ((j < k && dist.get(k, j) == UCHAR_MAX) || (j > k && dist.get(j, k) == UCHAR_MAX)))
            continue;
        unsigned char d = static_cast<unsigned char>(get_dist(i, k) * UCHAR_MAX);
        dist.set(i, k, d);
        if (d > 0) pq_.push({d, i, k});
    }
    for (int k = i + 1; k < ntrans; k++) {
        if (dist.no_link(k, i) && dist.no_link(k, j))
            continue;
        if (dist.get(k, i) == UCHAR_MAX && dist.get(k, j) == UCHAR_MAX)
            continue;
        unsigned char d = static_cast<unsigned char>(get_dist(i, k) * UCHAR_MAX);
        dist.set(k, i, d);
        if (d > 0) pq_.push({d, k, i});
    }

    // Remove j's entries
    for (int k = 0; k < j; k++) dist.remove(j, k);
    for (int k = j + 1; k < ntrans; k++) dist.remove(k, j);
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

void Cluster::cluster(map<float, string> &thresholds) {
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
            output_clusters(itr_d->second);
            itr_d++;
        }
        if (itr_d == thresholds.end() || distance == 1.0f) break;

        merge(i, j);
    }

    // Report final clustering for any remaining thresholds
    while (itr_d != thresholds.end()) {
        output_clusters(itr_d->second);
        itr_d++;
    }
}


// ── Output cluster assignments and counts ───────────────────────────

void Cluster::output_clusters(const string &threshold) {
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

    // Write atomically under critical section using persistent file handles
    #pragma omp critical(file_io)
    {
        static string last_counts_fn, last_cluster_fn;
        static ofstream countsFile, clusterFile;

        string counts_fn  = file_prefix + string(file_counts) + threshold + string(file_ext);
        string cluster_fn = file_prefix + string(file_clusters) + threshold + string(file_ext);

        if (counts_fn != last_counts_fn) {
            if (countsFile.is_open()) countsFile.close();
            countsFile.open(counts_fn, ios_base::app);
            last_counts_fn = counts_fn;
        }
        if (cluster_fn != last_cluster_fn) {
            if (clusterFile.is_open()) clusterFile.close();
            clusterFile.open(cluster_fn, ios_base::app);
            last_cluster_fn = cluster_fn;
        }

        countsFile  << counts_buf.str();  countsFile.flush();
        clusterFile << cluster_buf.str(); clusterFile.flush();
    }
}


// ── Initialise the distance matrix ──────────────────────────────────

void Cluster::initialise_matrix() {
    const int ntrans   = n_trans();
    const int nsamples = Transcript::samples;

    // Assign positional indices to transcripts
    for (int i = 0; i < ntrans; i++)
        get_tran(i)->pos(i);

    dist.set_size(ntrans);

    read_groups.resize(nsamples);
    for (int s = 0; s < nsamples; s++)
        read_groups[s].resize(ntrans);
    read_group_sizes.resize(ntrans);
    for (int t = 0; t < ntrans; t++)
        read_group_sizes[t].resize(nsamples, 0);

    // Build read-group lists and collect non-zero distance pairs
    vector<std::pair<int, int>> nz_pairs;
    for (int r = 0; r < n_reads(); r++) {
        Read *read   = get_read(r);
        int   sample = read->get_sample();

        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            for (auto t2 = read->align_begin(); t2 != t1; ++t2) {
                int j = (*t2)->pos();
                nz_pairs.emplace_back(std::max(i, j), std::min(i, j));
            }
            read_groups[sample][i].push_back(r);
            read_group_sizes[i][sample] += read->get_weight();
        }
    }

    // Each transcript gets its own initial cluster group
    groups.resize(ntrans);
    for (int n = 0; n < ntrans; n++)
        groups[n].push_back(n);

    // Deduplicate pairs to avoid redundant distance calculations
    std::ranges::sort(nz_pairs);
    auto [first, last] = std::ranges::unique(nz_pairs);
    nz_pairs.erase(first, last);

    // Compute distances and seed the max-heap
    dist.reserve(nz_pairs.size());
    for (auto &p : nz_pairs) {
        unsigned char d = static_cast<unsigned char>(get_dist(p.first, p.second) * UCHAR_MAX);
        dist.set(p.first, p.second, d);
        if (d > 0) pq_.push({d, p.first, p.second});
    }
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
