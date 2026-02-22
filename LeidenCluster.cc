// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
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
#include <string>
#include <algorithm>
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


// ── Build igraph from shared-read overlaps ──────────────────────────
// Reuses Cluster::get_dist() which already includes the LRT
// expression test (differential pairs get distance 0 → excluded).

static void build_graph(Cluster *c,
                        igraph_t *graph,
                        igraph_vector_t *weights,
                        int ntrans) {
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
    vector<int> edge_src, edge_dst;
    vector<double> edge_wt;
    edge_src.reserve(nz_pairs.size());
    edge_dst.reserve(nz_pairs.size());
    edge_wt.reserve(nz_pairs.size());

    for (auto &[hi, lo] : nz_pairs) {
        float d = c->get_dist(hi, lo);
        if (d > 0.0f) {
            edge_src.push_back(hi);
            edge_dst.push_back(lo);
            edge_wt.push_back(static_cast<double>(d));
        }
    }

    // Build igraph edge vector (flat: src0, dst0, src1, dst1, ...)
    igraph_integer_t nedges = static_cast<igraph_integer_t>(edge_src.size());
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 2 * nedges);
    for (igraph_integer_t e = 0; e < nedges; e++) {
        VECTOR(edges)[2 * e]     = edge_src[e];
        VECTOR(edges)[2 * e + 1] = edge_dst[e];
    }

    igraph_create(graph, &edges, ntrans, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges);

    // Copy edge weights
    igraph_vector_init(weights, nedges);
    for (igraph_integer_t e = 0; e < nedges; e++)
        VECTOR(*weights)[e] = edge_wt[e];
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

    // Build inverted index: read → list of communities it appears in
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
                          const vector<int> &membership) {
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
        string counts_fn  = Cluster::file_prefix + string(Cluster::file_counts)
                          + threshold + string(Cluster::file_ext);
        string cluster_fn = Cluster::file_prefix + string(Cluster::file_clusters)
                          + threshold + string(Cluster::file_ext);

        ofstream countsFile(counts_fn, ios_base::app);
        countsFile << counts_buf.str();

        ofstream clusterFile(cluster_fn, ios_base::app);
        clusterFile << cluster_buf.str();
    }
}


// ── Main Leiden entry point ─────────────────────────────────────────

void cluster_leiden(Cluster *c,
                    map<float, string> &thresholds) {
    const int ntrans = c->n_trans();

    // Always set up read-group data (needed for get_dist and count output)
    c->setup_read_groups();

    if (ntrans <= 1) {
        // Single transcript: just output it directly for each threshold
        for (auto &[thr, label] : thresholds) {
            vector<int> membership(ntrans, 0);
            output_leiden(c, label, 1, membership);
        }
        return;
    }

    if (ntrans > 1000) {
        #pragma omp critical(print)
        cout << "Leiden: cluster with " << ntrans
             << " transcripts" << endl;
    }

    // Seed per-thread RNG to avoid races on the global default RNG.
    // igraph_rng_default() returns thread-local storage when IGRAPH_THREAD_SAFE=1.
    igraph_rng_seed(igraph_rng_default(),
                    static_cast<igraph_uint_t>(c->get_id() + 42));

    // Build the shared-read graph (LRT pre-filtered)
    igraph_t graph;
    igraph_vector_t weights;
    build_graph(c, &graph, &weights, ntrans);

    // Run Leiden at each threshold's resolution
    for (auto &[thr, label] : thresholds) {
        igraph_vector_int_t igraph_membership;
        igraph_vector_int_init(&igraph_membership, ntrans);
        igraph_int_t nb_clusters = 0;
        igraph_real_t quality = 0;

        // CPM resolution = threshold value (edge weights are in [0,1])
        double resolution = static_cast<double>(thr);

        igraph_community_leiden_simple(
            &graph,
            &weights,
            IGRAPH_LEIDEN_OBJECTIVE_CPM,
            resolution,
            /* beta= */ 0.01,
            /* start= */ 0,          // fresh partition
            /* n_iterations= */ -1,  // iterate until convergence
            &igraph_membership,
            &nb_clusters,
            &quality
        );

        // Convert igraph membership to std::vector
        int ncomm = static_cast<int>(nb_clusters);
        vector<int> membership(ntrans);
        for (int i = 0; i < ntrans; i++)
            membership[i] = static_cast<int>(VECTOR(igraph_membership)[i]);

        igraph_vector_int_destroy(&igraph_membership);

        output_leiden(c, label, ncomm, membership);
    }

    // Cleanup
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
}

#endif // HAVE_IGRAPH
