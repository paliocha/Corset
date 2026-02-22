// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
//
// Features:
//   - Soft LRT weighting: sigmoid decay around D_cut (--lrt-softness)
//   - Adaptive kNN sparsification: connectivity-constrained minimum-k
//     per super-cluster (--knn auto), or fixed global k (--knn <int>)
//
// Author: Martin Paliocha, martin.paliocha@nmbu.no
// Created 22 February 2026

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
#include <queue>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ranges>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::ios_base;


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


// ── Collect edges from shared reads ─────────────────────────────────
// Reuses Cluster::get_dist() which includes the LRT expression test.
// Passes lrt_softness through so Leiden can use sigmoid decay.

static EdgeList collect_edges(Cluster *c, int ntrans) {
    // Collect non-zero transcript pairs from shared reads
    vector<std::pair<int, int>> nz_pairs;
    for (int r = 0; r < c->n_reads(); r++) {
        Read *read = c->get_read(r);
        for (auto t1 = read->align_begin(); t1 != read->align_end(); ++t1) {
            int i = (*t1)->pos();
            for (auto t2 = read->align_begin(); t2 != t1; ++t2) {
                int j = (*t2)->pos();
                nz_pairs.emplace_back(std::max(i, j), std::min(i, j));
            }
        }
    }

    // Deduplicate
    std::ranges::sort(nz_pairs);
    auto [first, last] = std::ranges::unique(nz_pairs);
    nz_pairs.erase(first, last);

    // Compute distances and collect edges with non-zero weight
    EdgeList edges;
    edges.reserve(nz_pairs.size());

    const float softness = Cluster::lrt_softness;
    for (auto &[hi, lo] : nz_pairs) {
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


// ── Count connected components via BFS ──────────────────────────────
// O(V + E), used by the connectivity-constrained k search.

static int count_components(const EdgeList &edges, int ntrans) {
    vector<vector<int>> adj(ntrans);
    for (size_t e = 0; e < edges.size(); e++) {
        adj[edges.src[e]].push_back(edges.dst[e]);
        adj[edges.dst[e]].push_back(edges.src[e]);
    }

    vector<bool> visited(ntrans, false);
    int components = 0;
    std::queue<int> q;
    for (int n = 0; n < ntrans; n++) {
        if (visited[n]) continue;
        components++;
        q.push(n);
        visited[n] = true;
        while (!q.empty()) {
            int cur = q.front();
            q.pop();
            for (int nb : adj[cur]) {
                if (!visited[nb]) {
                    visited[nb] = true;
                    q.push(nb);
                }
            }
        }
    }
    return components;
}


// ── Find minimum k preserving graph connectivity ────────────────────
// Binary search in [1, sqrt(n)] for the smallest k where the number
// of connected components does not exceed target_components.

static int find_optimal_k(const EdgeList &edges, int ntrans,
                          int target_components) {
    int lo = 1;
    int hi = static_cast<int>(std::round(std::sqrt(ntrans)));
    int best_k = hi;

    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        EdgeList filtered = apply_knn(edges, ntrans, mid);
        int comp = count_components(filtered, ntrans);

        if (comp <= target_components) {
            best_k = mid;
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    return best_k;
}


// ── Run Leiden and return membership + community count ───────────────

struct LeidenResult {
    vector<int> membership;
    int         ncommunities;
    double      quality;
};

static LeidenResult run_leiden(igraph_t *graph, igraph_vector_t *weights,
                               double resolution, int ntrans) {
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
        /* n_iterations= */ -1,
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

    // Randomly assign each read to one of its communities
    vector<int> counts(ncommunities, 0);
    for (int rd = 0; rd < nreads; rd++) {
        int n_align = static_cast<int>(read_to_comm[rd].size());
        if (n_align == 0) continue;
        int weight = c->get_read(rd)->get_weight();
        for (; weight > 0; weight--)
            counts[read_to_comm[rd][rand_r(&seed) % n_align]]++;
    }
    return counts;
}


// ── Output Leiden communities in Corset format ──────────────────────

static void output_leiden(Cluster *c,
                          const string &threshold,
                          int ncommunities,
                          const vector<int> &membership,
                          const string &method_tag) {
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

    // Write atomically under critical section
    #pragma omp critical(file_io)
    {
        string counts_fn  = Cluster::file_prefix + method_tag
                          + string(Cluster::file_counts)
                          + threshold + string(Cluster::file_ext);
        string cluster_fn = Cluster::file_prefix + method_tag
                          + string(Cluster::file_clusters)
                          + threshold + string(Cluster::file_ext);

        ofstream countsFile(counts_fn, ios_base::app);
        countsFile << counts_buf.str();

        ofstream clusterFile(cluster_fn, ios_base::app);
        clusterFile << cluster_buf.str();
    }
}


// ── Main Leiden entry point ─────────────────────────────────────────

void cluster_leiden(Cluster *c,
                    map<float, string> &thresholds,
                    const string &method_tag) {
    const int ntrans = c->n_trans();

    // Always set up read-group data (needed for get_dist and count output)
    c->setup_read_groups();

    if (ntrans <= 1) {
        for (auto &[thr, label] : thresholds) {
            vector<int> membership(ntrans, 0);
            output_leiden(c, label, 1, membership, method_tag);
        }
        return;
    }

    if (ntrans > 1000) {
        #pragma omp critical(print)
        cout << "Leiden: super-cluster " << c->get_id()
             << " with " << ntrans << " transcripts" << endl;
    }

    // Seed per-thread RNG for deterministic Leiden results
    igraph_rng_seed(igraph_rng_default(),
                    static_cast<igraph_uint_t>(c->get_id() + 42));

    // Collect all edges (expensive: calls get_dist for every shared-read pair)
    EdgeList edges = collect_edges(c, ntrans);

    // ── kNN sparsification ──────────────────────────────────────────
    const int knn = Cluster::knn;

    if (knn >= 0 && ntrans > 2) {
        int effective_k;

        if (knn == 0) {
            // Auto mode: connectivity-constrained minimum k
            // Step 1: run Leiden on full graph to get target community count
            igraph_t full_graph;
            igraph_vector_t full_weights;
            build_igraph(edges, ntrans, &full_graph, &full_weights);

            double first_resolution = static_cast<double>(thresholds.begin()->first);
            LeidenResult full_result = run_leiden(&full_graph, &full_weights,
                                                  first_resolution, ntrans);
            int target_components = full_result.ncommunities;

            igraph_vector_destroy(&full_weights);
            igraph_destroy(&full_graph);

            // Step 2: binary search for minimum k preserving connectivity
            effective_k = find_optimal_k(edges, ntrans, target_components);

            if (ntrans > 1000) {
                size_t full_edges = edges.size();
                EdgeList filtered = apply_knn(edges, ntrans, effective_k);
                #pragma omp critical(print)
                cout << "  kNN auto: k=" << effective_k
                     << " (target " << target_components
                     << " communities, " << full_edges
                     << " -> " << filtered.size() << " edges)" << endl;
                edges = std::move(filtered);
            } else {
                edges = apply_knn(edges, ntrans, effective_k);
            }
        } else {
            // Fixed k mode
            effective_k = std::min(knn, ntrans - 1);
            edges = apply_knn(edges, ntrans, effective_k);
        }
    }

    // Build igraph from (potentially kNN-filtered) edges
    igraph_t graph;
    igraph_vector_t weights;
    build_igraph(edges, ntrans, &graph, &weights);

    // Run Leiden at each threshold's resolution
    for (auto &[thr, label] : thresholds) {
        double resolution = static_cast<double>(thr);
        LeidenResult result = run_leiden(&graph, &weights, resolution, ntrans);
        output_leiden(c, label, result.ncommunities, result.membership, method_tag);
    }

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
}

#endif // HAVE_IGRAPH
