// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
//
// Features:
//   - Soft LRT weighting: sigmoid decay around D_cut (--lrt-softness)
//   - Adaptive kNN sparsification: quality-preserving minimum-k per
//     super-cluster (--knn auto), or fixed global k (--knn <int>)
//   - Two-phase edge construction: hash-based pair enumeration with
//     cheap proxy ranking, get_dist only on kNN candidates
//
// Author: Martin Paliocha, martin.paliocha@nmbu.no
// Created 22 February 2026, last modified 23 February 2026

#ifdef HAVE_IGRAPH

#include <LeidenCluster.h>
#include <Cluster.h>
#include <Transcript.h>
#include <Read.h>

#include <igraph/igraph.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ranges>
#include <omp.h>
#include <FlatDistMap.h>
#include <Progress.h>

using std::cout;
using std::endl;
using std::ios_base;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::vector;


// ── Edge list: raw edges before igraph construction ─────────────────
// Separated from igraph creation so kNN filtering can operate on the
// flat vectors before building the graph object.

struct EdgeList {
    vector<int>    src;
    vector<int>    dst;
    vector<double> wt;

    void reserve(size_t n) { src.reserve(n); dst.reserve(n); wt.reserve(n); }
    [[nodiscard]] size_t size() const { return src.size(); }
};


// ── Direct pairwise edge collection ─────────────────────────────────
// Computes get_dist() for all C(ntrans, 2) pairs directly.  Avoids the
// O(nreads × alignments²) read-based pair enumeration that dominates
// for small super-clusters with highly multi-mapped reads (common in
// de novo transcriptomes).  Cost: O(ntrans² × nsamples × avg_reads)
// for the sorted-list intersections inside get_dist().
//
// Used when the cost of direct pairwise (C(ntrans,2) get_dist calls)
// is cheaper than read-based enumeration (Σ C(alignments,2) hash ops).

static constexpr int DIRECT_EDGE_MAX = 20000;

// Estimate cost of read-based pair enumeration: Σ C(alignments_r, 2).
// Returns true when direct pairwise (C(ntrans,2) get_dist calls) is cheaper.
static bool prefer_direct(Cluster *c, int ntrans) {
    if (ntrans > DIRECT_EDGE_MAX) return false;
    int64_t direct_cost = static_cast<int64_t>(ntrans) * (ntrans - 1) / 2;
    int64_t read_cost = 0;
    for (int r = 0; r < c->n_reads(); r++) {
        int a = c->get_read(r)->alignments();
        read_cost += static_cast<int64_t>(a) * (a - 1) / 2;
        if (read_cost > direct_cost) return true;  // early out
    }
    return false;
}

static EdgeList collect_edges_direct(Cluster *c, int ntrans) {
    const float softness = Cluster::lrt_softness;
    EdgeList edges;
    for (int i = 0; i < ntrans; i++) {
        for (int j = 0; j < i; j++) {
            float d = c->get_dist(i, j, softness);
            if (d > 0.0f) {
                edges.src.push_back(i);
                edges.dst.push_back(j);
                edges.wt.push_back(static_cast<double>(d));
            }
        }
    }
    return edges;
}


// ── Collect edges from shared reads (hash-based dedup) ──────────────
// Uses FlatDistMap as a seen-set for O(1) dedup per raw pair.
// Replaces the old sort+unique approach: O(N) vs O(N log N), and
// peak memory is O(U) not O(N) where U << N = unique vs raw pairs.

static EdgeList collect_edges(Cluster *c, int ntrans) {
    const auto nt64 = static_cast<uint64_t>(ntrans);
    const float softness = Cluster::lrt_softness;

    FlatDistMap seen;
    seen.reserve(ntrans * 4);
    EdgeList edges;

    for (int r = 0; r < c->n_reads(); r++) {
        Read *read = c->get_read(r);
        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            for (auto t2 = read->align_begin(); t2 != t1; ++t2) {
                int j = (*t2)->pos();
                int hi = std::max(i, j), lo = std::min(i, j);
                uint64_t key = nt64 * hi + lo;
                if (!seen.contains(key)) {
                    seen.set(key, 1);
                    float d = c->get_dist(hi, lo, softness);
                    if (d > 0.0f) {
                        edges.src.push_back(hi);
                        edges.dst.push_back(lo);
                        edges.wt.push_back(static_cast<double>(d));
                    }
                }
            }
        }
    }
    return edges;
}


// ── Two-phase edge collection for kNN modes ─────────────────────────
// Phase 1: Enumerate unique pairs via hash, count occurrences as a
//          cheap proxy for distance (more shared reads → likely higher
//          distance).  Uses FlatDistMap saturating counter (0-255).
// Phase 1b: Per-node, keep top-k_max neighbors by proxy score.
//           Symmetric kNN: keep if EITHER endpoint selects it.
// Phase 2: Compute get_dist() only on selected candidates.
// Avoids calling the expensive get_dist() on pairs that kNN will
// discard — decisive for dense super-clusters.

static EdgeList collect_edges_knn(Cluster *c, int ntrans, int k_max) {
    const auto nt64 = static_cast<uint64_t>(ntrans);
    const float softness = Cluster::lrt_softness;

    // ── Phase 1: Hash-based pair enumeration with occurrence counting ──
    FlatDistMap pair_counts;
    pair_counts.reserve(ntrans * 4);
    vector<vector<int>> node_adj(ntrans);

    for (int r = 0; r < c->n_reads(); r++) {
        Read *read = c->get_read(r);
        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            for (auto t2 = read->align_begin(); t2 != t1; ++t2) {
                int j = (*t2)->pos();
                int hi = std::max(i, j), lo = std::min(i, j);
                uint64_t key = nt64 * hi + lo;
                unsigned char cur = pair_counts.get(key);
                if (cur == 0) {
                    node_adj[hi].push_back(lo);
                    node_adj[lo].push_back(hi);
                }
                if (cur < 255) pair_counts.set(key, cur + 1);
            }
        }
    }

    // ── Phase 1b: Per-node top-k_max selection by proxy score ──
    // Build a deduped candidate list for Phase 2 so node_adj can be freed
    // eagerly — peak memory drops as adjacency shrinks and candidates grow.
    FlatDistMap seen_pair;
    seen_pair.reserve(static_cast<size_t>(ntrans) * k_max / 2);
    vector<std::pair<int, int>> candidates;
    candidates.reserve(static_cast<size_t>(ntrans) * k_max / 2);

    vector<std::pair<unsigned char, int>> scored;  // reused across nodes
    for (int n = 0; n < ntrans; n++) {
        auto &adj = node_adj[n];

        // Lambda: nominate a neighbor, dedup via seen_pair
        auto nominate = [&](int nb) {
            int hi = std::max(n, nb), lo = std::min(n, nb);
            uint64_t key = nt64 * hi + lo;
            if (!seen_pair.contains(key)) {
                seen_pair.set(key, 1);
                candidates.emplace_back(hi, lo);
            }
        };

        if (static_cast<int>(adj.size()) <= k_max) {
            for (int nb : adj) nominate(nb);
        } else {
            scored.clear();
            for (int nb : adj) {
                uint64_t key = nt64 * std::max(n, nb) + std::min(n, nb);
                scored.emplace_back(pair_counts.get(key), nb);
            }
            std::nth_element(scored.begin(), scored.begin() + k_max, scored.end(),
                [](auto &a, auto &b) { return a.first > b.first; });
            for (int i = 0; i < k_max; i++)
                nominate(scored[i].second);
        }

        adj = {};  // free adjacency eagerly
    }

    // Release Phase 1 structures (no longer needed)
    node_adj = {};
    pair_counts = FlatDistMap{};
    seen_pair = FlatDistMap{};

    // ── Phase 2: Compute get_dist() on selected candidates only ──
    EdgeList edges;
    edges.reserve(candidates.size());

    for (auto [hi, lo] : candidates) {
        float d = c->get_dist(hi, lo, softness);
        if (d > 0.0f) {
            edges.src.push_back(hi);
            edges.dst.push_back(lo);
            edges.wt.push_back(static_cast<double>(d));
        }
    }

    return edges;
}


// ── Build igraph from edge list ─────────────────────────────────────

static void build_igraph(const EdgeList &edges, int ntrans,
                         igraph_t *graph, igraph_vector_t *weights) {
    igraph_integer_t nedges = static_cast<igraph_integer_t>(edges.size());
    igraph_vector_int_t igraph_edges;
    igraph_vector_int_init(&igraph_edges, 2 * nedges);
    for (igraph_integer_t e = 0; e < nedges; e++) {
        VECTOR(igraph_edges)[2 * e]     = edges.src[e];
        VECTOR(igraph_edges)[2 * e + 1] = edges.dst[e];
    }

    igraph_create(graph, &igraph_edges, ntrans, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&igraph_edges);

    igraph_vector_init(weights, nedges);
    for (igraph_integer_t e = 0; e < nedges; e++)
        VECTOR(*weights)[e] = edges.wt[e];
}


// ── kNN sparsification ──────────────────────────────────────────────
// Keep only the k highest-weight edges per node (symmetric: keep if
// EITHER endpoint claims the edge in its top-k).

static EdgeList apply_knn(const EdgeList &in, int ntrans, int k) {
    const size_t nedges = in.size();

    // Per-node: collect (weight, edge_index) pairs
    vector<vector<std::pair<double, size_t>>> node_edges(ntrans);
    for (size_t e = 0; e < nedges; e++) {
        node_edges[in.src[e]].emplace_back(in.wt[e], e);
        node_edges[in.dst[e]].emplace_back(in.wt[e], e);
    }

    // For each node, mark top-k edges by weight
    vector<bool> keep(nedges, false);
    for (int n = 0; n < ntrans; n++) {
        auto &ne = node_edges[n];
        if (static_cast<int>(ne.size()) <= k) {
            for (auto &[w, idx] : ne) keep[idx] = true;
            continue;
        }
        // Partial sort: top-k by weight (descending)
        std::nth_element(ne.begin(), ne.begin() + k, ne.end(),
                         [](auto &a, auto &b) { return a.first > b.first; });
        for (int i = 0; i < k; i++)
            keep[ne[i].second] = true;
    }

    // Compact kept edges
    EdgeList out;
    out.reserve(nedges);
    for (size_t e = 0; e < nedges; e++) {
        if (keep[e]) {
            out.src.push_back(in.src[e]);
            out.dst.push_back(in.dst[e]);
            out.wt.push_back(in.wt[e]);
        }
    }
    return out;
}


// ── Leiden result container ─────────────────────────────────────────

struct LeidenResult {
    vector<int> membership;
    int         ncommunities;
    double      quality;
};


// ── Compare two partitions via Normalized Mutual Information ────────
// Wraps igraph_compare_communities().  Returns NMI in [0, 1] where
// 1 = identical partitions.

static double compute_nmi(const vector<int> &a, const vector<int> &b, int n) {
    igraph_vector_int_t va, vb;
    igraph_vector_int_init(&va, n);
    igraph_vector_int_init(&vb, n);
    for (int i = 0; i < n; i++) {
        VECTOR(va)[i] = a[i];
        VECTOR(vb)[i] = b[i];
    }
    igraph_real_t nmi = 0;
    igraph_compare_communities(&va, &vb, &nmi, IGRAPH_COMMCMP_NMI);
    igraph_vector_int_destroy(&vb);
    igraph_vector_int_destroy(&va);
    return static_cast<double>(nmi);
}


// ── Quality-preserving kNN search ───────────────────────────────────
// Binary search in [floor, k_max] for the minimum k where:
//   (1) CPM quality >= QUALITY_THRESHOLD * Q_ref
//   (2) NMI(partition_knn, partition_ref) >= NMI_THRESHOLD
// Uses fast Leiden (2 iterations) during search for efficiency.
// The k floor prevents the degenerate k=1 case that the old
// connectivity-based search suffered from.

static constexpr double KNN_QUALITY_THRESHOLD = 0.8;
static constexpr double KNN_NMI_THRESHOLD     = 0.7;

// Forward declaration (defined below).
static LeidenResult run_leiden(igraph_t *graph, igraph_vector_t *weights,
                               double resolution, int ntrans,
                               int n_iterations);

// sc_id >= 0: log each binary search probe (for large SCs).
// sc_id < 0:  quiet mode.
static int find_optimal_k(const EdgeList &edges, int ntrans,
                          double resolution,
                          const LeidenResult &ref_result,
                          int sc_id = -1) {
    int k_floor = std::max(3, static_cast<int>(std::ceil(std::log2(ntrans + 1))));
    int k_max   = static_cast<int>(std::ceil(std::sqrt(ntrans)));
    int lo      = k_floor;
    int hi      = k_max;
    int best_k  = k_max;

    double Q_ref = ref_result.quality;

    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        EdgeList filtered = apply_knn(edges, ntrans, mid);

        igraph_t graph;
        igraph_vector_t weights;
        build_igraph(filtered, ntrans, &graph, &weights);

        // Fast Leiden: 2 iterations (not full convergence)
        LeidenResult knn_result = run_leiden(&graph, &weights,
                                             resolution, ntrans, 2);

        igraph_vector_destroy(&weights);
        igraph_destroy(&graph);

        // Quality check: handles positive, negative, and near-zero Q_ref.
        // Multiplying threshold by Q_ref (rather than dividing quality)
        // preserves the inequality direction when Q_ref is negative.
        bool quality_ok = (std::abs(Q_ref) > 1e-15)
            ? knn_result.quality >= KNN_QUALITY_THRESHOLD * Q_ref
            : true;

        // NMI structural similarity check
        double nmi = compute_nmi(ref_result.membership,
                                 knn_result.membership, ntrans);

        bool ok = quality_ok && (nmi >= KNN_NMI_THRESHOLD);

        if (sc_id >= 0) {
            double q_ratio = (std::abs(Q_ref) > 1e-15)
                ? knn_result.quality / Q_ref : 1.0;
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                          "      probe k=%d: Q=%.2f*Q_ref, NMI=%.2f %s",
                          mid, q_ratio, nmi,
                          ok ? "\xe2\x9c\x93" : "\xe2\x9c\x97");
            #pragma omp critical(print)
            cout << progress::ansi::dim(buf) << endl;
        }

        if (ok) {
            best_k = mid;
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }

    return best_k;
}


// ── Run Leiden and return membership + community count ───────────────
// n_iterations = -1 for full convergence, or a small positive value
// (e.g. 2) for a fast approximation during kNN search.

static LeidenResult run_leiden(igraph_t *graph, igraph_vector_t *weights,
                               double resolution, int ntrans,
                               int n_iterations = -1) {
    igraph_vector_int_t igraph_membership;
    igraph_vector_int_init(&igraph_membership, ntrans);
    igraph_int_t nb_clusters = 0;
    igraph_real_t quality = 0;

    igraph_community_leiden_simple(
        graph, weights,
        IGRAPH_LEIDEN_OBJECTIVE_CPM,
        resolution,
        /* beta= */ 0.01,
        /* start= */ 0,
        static_cast<igraph_int_t>(n_iterations),
        &igraph_membership,
        &nb_clusters,
        &quality
    );

    LeidenResult result;
    result.ncommunities = static_cast<int>(nb_clusters);
    result.quality      = static_cast<double>(quality);
    result.membership.resize(ntrans);
    for (int i = 0; i < ntrans; i++)
        result.membership[i] = static_cast<int>(VECTOR(igraph_membership)[i]);

    igraph_vector_int_destroy(&igraph_membership);
    return result;
}


// ── Count reads per Leiden community for a given sample ─────────────
// Mirrors Cluster::get_counts() but uses Leiden community assignments
// instead of hierarchical merge groups.

static vector<int> get_counts_leiden(Cluster *c,
                                     int sample,
                                     int ncommunities,
                                     const vector<int> &membership) {
    const int nreads  = c->n_reads();
    const int ntrans  = c->n_trans();
    const auto &rg    = c->get_read_groups();

    // Thread-safe RNG seed
    unsigned int seed = static_cast<unsigned int>(c->get_id() * 17 + sample * 31 + 42);

    // Build inverted index: read -> list of communities it appears in
    vector<vector<int>> read_to_comm(nreads);
    for (int t = 0; t < ntrans; t++) {
        int comm = membership[t];
        for (int idx : rg[sample][t])
            read_to_comm[idx].push_back(comm);
    }

    // Deduplicate
    for (int rd = 0; rd < nreads; rd++) {
        auto &v = read_to_comm[rd];
        if (v.size() > 1) {
            std::ranges::sort(v);
            auto [f, l] = std::ranges::unique(v);
            v.erase(f, l);
        }
    }

    // Assign each read's weight to one of its communities.
    // When a read maps to exactly one community (common case), add weight
    // directly — avoids an O(weight) loop per equivalence class.
    vector<int> counts(ncommunities, 0);
    for (int rd = 0; rd < nreads; rd++) {
        int n_align = static_cast<int>(read_to_comm[rd].size());
        if (n_align == 0) continue;
        int weight = c->get_read(rd)->get_weight();
        if (n_align == 1) {
            counts[read_to_comm[rd][0]] += weight;
        } else {
            // O(n_align) even-division instead of O(weight) random loop.
            // Expected value identical; for large weights the stochastic
            // variance is negligible (< sqrt(weight)).
            int per = weight / n_align;
            int remainder = weight % n_align;
            for (int c = 0; c < n_align; c++)
                counts[read_to_comm[rd][c]] += per;
            for (int r = 0; r < remainder; r++)
                counts[read_to_comm[rd][rand_r(&seed) % n_align]]++;
        }
    }

    // Add direct counts per Leiden community
    for (int t = 0; t < ntrans; t++) {
        Transcript *tr = c->get_tran(t);
        if (tr->has_direct_counts())
            counts[membership[t]] += tr->get_direct_count(sample);
    }

    return counts;
}


// ── Output Leiden communities in Corset format ──────────────────────

static void output_leiden(Cluster *c,
                          const string &threshold,
                          int ncommunities,
                          const vector<int> &membership,
                          const string &method_tag,
                          map<string, string> *file_buf = nullptr) {
    const int ntrans = c->n_trans();

    // Compute counts per community per sample
    vector<vector<int>> counts;
    counts.reserve(Transcript::samples);
    for (int s = 0; s < Transcript::samples; s++)
        counts.push_back(get_counts_leiden(c, s, ncommunities, membership));

    // Precompute which communities have at least one member (O(ntrans))
    vector<bool> has_members(ncommunities, false);
    for (int t = 0; t < ntrans; t++)
        has_members[membership[t]] = true;

    // Buffer counts output
    ostringstream counts_buf;
    for (int comm = 0; comm < ncommunities; comm++) {
        if (!has_members[comm]) continue;

        counts_buf << Cluster::cluster_id_prefix << c->get_id()
                   << Cluster::cluster_id_joiner << comm;
        for (int s = 0; s < Transcript::samples; s++)
            counts_buf << "\t" << counts[s][comm];
        counts_buf << "\n";
    }

    // Buffer cluster assignments
    ostringstream cluster_buf;
    for (int t = 0; t < ntrans; t++) {
        cluster_buf << c->get_tran(t)->get_name()
                    << "\t" << Cluster::cluster_id_prefix << c->get_id()
                    << Cluster::cluster_id_joiner << membership[t] << "\n";
    }

    string counts_fn  = Cluster::file_prefix + method_tag
                      + string(Cluster::file_counts)
                      + threshold + string(Cluster::file_ext);
    string cluster_fn = Cluster::file_prefix + method_tag
                      + string(Cluster::file_clusters)
                      + threshold + string(Cluster::file_ext);

    if (file_buf) {
        (*file_buf)[counts_fn] += counts_buf.str();
        (*file_buf)[cluster_fn] += cluster_buf.str();
    } else {
        #pragma omp critical(file_io)
        {
            ofstream countsFile(counts_fn, ios_base::app);
            countsFile << counts_buf.str();

            ofstream clusterFile(cluster_fn, ios_base::app);
            clusterFile << cluster_buf.str();
        }
    }
}


// ── Main Leiden entry point ─────────────────────────────────────────

int cluster_leiden(Cluster *c,
                   map<float, string> &thresholds,
                   const string &method_tag,
                   vector<int> *cluster_sizes,
                   map<string, string> *file_buf) {
    const int ntrans = c->n_trans();
    const int nreads = c->n_reads();

    // Always set up read-group data (needed for get_dist and count output)
    double t_setup = omp_get_wtime();
    c->setup_read_groups();
    double t_after_setup = omp_get_wtime();

    if (ntrans <= 1) {
        for (auto &[thr, label] : thresholds) {
            vector<int> membership(ntrans, 0);
            output_leiden(c, label, 1, membership, method_tag, file_buf);
        }
        if (cluster_sizes && ntrans > 0)
            cluster_sizes->assign(1, ntrans);
        return 0;
    }

    const bool large = ntrans > 1000;
    if (large) {
        #pragma omp critical(print)
        cout << progress::ansi::cyan(
                "  \xe2\x96\xb8 SC " + std::to_string(c->get_id()) + ": "
                + progress::format_count(ntrans) + " transcripts (Leiden)")
             << endl;
    }

    double t0 = t_after_setup;
    const int knn = Cluster::knn;
    const double first_resolution = static_cast<double>(thresholds.begin()->first);
    EdgeList edges;
    double t_edges = 0, t_ref = 0, t_search = 0;

    // Decide edge-collection strategy: direct pairwise (O(ntrans²) get_dist
    // calls) vs read-based enumeration (O(Σ C(alignments,2)) hash ops).
    // Direct wins when reads are heavily multi-mapped (de novo transcriptomes).
    const bool use_direct = ntrans > 2 && prefer_direct(c, ntrans);

    if (knn == 0 && use_direct) {
        // ── kNN auto, direct pairwise cheaper than read-based enumeration ──
        // Full graph is tractable — skip kNN sparsification entirely.
        edges = collect_edges_direct(c, ntrans);
        t_edges = omp_get_wtime() - t0;

    } else if (knn == 0 && ntrans > 2) {
        // ── kNN auto: two-phase construction + quality-preserving search ──
        int k_max = static_cast<int>(std::ceil(std::sqrt(ntrans)));

        // Two-phase: cheap proxy ranking then get_dist on K_max candidates
        edges = collect_edges_knn(c, ntrans, k_max);
        double t1 = omp_get_wtime();
        t_edges = t1 - t0;

        // kNN search only meaningful when edges exist (all-zero get_dist
        // produces an empty graph → singletons regardless of k)
        if (edges.size() > 0) {
            // Reference Leiden run on K_max-sparse graph
            igraph_t ref_graph;
            igraph_vector_t ref_weights;
            build_igraph(edges, ntrans, &ref_graph, &ref_weights);
            LeidenResult ref_result = run_leiden(&ref_graph, &ref_weights,
                                                 first_resolution, ntrans);
            igraph_vector_destroy(&ref_weights);
            igraph_destroy(&ref_graph);
            double t2 = omp_get_wtime();
            t_ref = t2 - t1;

            // Quality-preserving binary search for minimum k
            int best_k = find_optimal_k(edges, ntrans, first_resolution,
                                        ref_result,
                                        large ? c->get_id() : -1);
            double t3 = omp_get_wtime();
            t_search = t3 - t2;

            size_t full_edges = edges.size();
            edges = apply_knn(edges, ntrans, best_k);

            if (large) {
                #pragma omp critical(print)
                {
                    cout << "    kNN auto: k=" << best_k
                         << " (Q_ref=" << ref_result.quality
                         << ", " << ref_result.ncommunities << " communities"
                         << ", " << progress::format_count(static_cast<int64_t>(full_edges))
                         << " -> " << progress::format_count(static_cast<int64_t>(edges.size()))
                         << " edges)" << endl;
                    cout << progress::ansi::dim(
                            "    [edges=" + progress::format_duration(t_edges)
                            + " ref=" + progress::format_duration(t_ref)
                            + " search=" + progress::format_duration(t_search) + "]")
                         << endl;
                }
            }
        }

    } else if (knn > 0 && use_direct) {
        // ── Fixed k, direct pairwise then kNN filter ──
        edges = collect_edges_direct(c, ntrans);
        edges = apply_knn(edges, ntrans, std::min(knn, ntrans - 1));
        t_edges = omp_get_wtime() - t0;

    } else if (knn > 0 && ntrans > 2) {
        // ── Fixed k: two-phase with generous K_max, then filter ──
        int k_max = std::min(2 * knn, ntrans - 1);
        edges = collect_edges_knn(c, ntrans, k_max);
        edges = apply_knn(edges, ntrans, std::min(knn, ntrans - 1));
        t_edges = omp_get_wtime() - t0;

    } else if (use_direct) {
        // ── No kNN, direct pairwise ──
        edges = collect_edges_direct(c, ntrans);
        t_edges = omp_get_wtime() - t0;

    } else {
        // ── No kNN: hash-dedup, all edges ──
        edges = collect_edges(c, ntrans);
        t_edges = omp_get_wtime() - t0;
    }

    int nedges = static_cast<int>(edges.size());

    // Build igraph from (potentially kNN-filtered) edges
    igraph_t graph;
    igraph_vector_t weights;
    build_igraph(edges, ntrans, &graph, &weights);

    // Run Leiden at each threshold's resolution
    double t_before_output = omp_get_wtime();
    for (auto &[thr, label] : thresholds) {
        double resolution = static_cast<double>(thr);
        LeidenResult result = run_leiden(&graph, &weights, resolution, ntrans);
        output_leiden(c, label, result.ncommunities, result.membership, method_tag, file_buf);

        // Capture cluster sizes from the last threshold
        if (cluster_sizes) {
            cluster_sizes->assign(result.ncommunities, 0);
            for (int m : result.membership)
                (*cluster_sizes)[m]++;
        }
    }
    double t_output = omp_get_wtime() - t_before_output;

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    double t_end = omp_get_wtime();
    double t_total = t_end - t_setup;

    if (large) {
        #pragma omp critical(print)
        cout << progress::ansi::green(
                "    SC " + std::to_string(c->get_id()) + " done: "
                + progress::format_duration(t_total) + " ("
                + progress::format_count(nedges) + " edges)")
             << endl;
    }

    // Diagnostic: log any SC taking >1s (catches slow small SCs)
    if (!large && t_total > 1.0) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "    \xe2\x9a\xa0 SC %d: %d transcripts, %d reads, %d edges, "
                      "%.1fs total [setup=%.1fs edges=%.1fs ref=%.1fs search=%.1fs output=%.1fs]",
                      c->get_id(), ntrans, nreads, nedges,
                      t_total, t_after_setup - t_setup, t_edges, t_ref, t_search, t_output);
        #pragma omp critical(print)
        cout << progress::ansi::yellow(buf) << endl;
    }

    return nedges;
}

#endif // HAVE_IGRAPH
