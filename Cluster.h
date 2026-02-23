// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// A Cluster (super-cluster) holds all transcripts that share
// at least one read.  Two clustering back-ends are available:
//
// 1. Hierarchical (default) — agglomerative clustering with a
//    log-likelihood ratio expression test.  Each merge step
//    walks per-group adjacency lists rather than scanning all
//    transcript groups (O(degree) per merge).  Closest pairs
//    are extracted via a lazy-deletion max-heap.  Distance
//    recomputation after merges is OpenMP-parallelised when
//    the merged group has many live neighbors.
//
// 2. Leiden (--algorithm leiden) — community detection via
//    igraph's CPM implementation.  See LeidenCluster.cc.
//
// Original author: Nadia Davidson
// Last modified 22 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <fstream>
#include <vector>
#include <queue>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <cstdint>
#include <climits>
#include <FlatDistMap.h>
#include <Read.h>

using group       = std::vector<std::vector<int>>;
using read_group  = std::vector<std::vector<std::vector<int>>>;

enum class ClusterAlgorithm { Hierarchical, Leiden, Both };

// Lazy-deletion max-heap entry for find_next_pair().
struct HeapEntry {
    unsigned char value;
    int i, j;   // invariant: i > j
    [[nodiscard]] constexpr bool operator<(const HeapEntry &o) const { return value < o.value; }
};

// Sparse distance matrix keyed by (i*ntrans + j), backed by a
// flat open-addressing hash map for ~2-3× better cache performance
// than std::unordered_map's chained buckets.
class DistanceMatrix {
    FlatDistMap dist_;
    uint64_t ntrans_ = 0;

public:
    void set_size(int n)    { ntrans_ = static_cast<uint64_t>(n); }
    void reserve(size_t n)  { dist_.reserve(n); }

    [[nodiscard]] unsigned char get(int i, int j) const {
        return dist_.get(ntrans_ * i + j);
    }
    void set(int i, int j, int value) {
        uint64_t key = ntrans_ * i + j;
        if (value != 0) dist_.set(key, static_cast<unsigned char>(value));
        else            dist_.erase(key);
    }
    [[nodiscard]] bool no_link(int i, int j) const { return !dist_.contains(ntrans_ * i + j); }
    void remove(int i, int j)                      { dist_.erase(ntrans_ * i + j); }
};


class Cluster {
    // --- Data ---
    std::vector<Transcript *> cluster_;
    std::vector<Read *>       read_;
    read_group read_groups;       // [sample][transcript][read indices]
    group      groups;            // [transcript] → list of grouped transcript indices
    group      read_group_sizes;  // [transcript][sample] → weighted read count

    DistanceMatrix dist;
    std::priority_queue<HeapEntry> pq_;

    // Adjacency lists: adj_[i] = sorted list of live neighbors with non-zero
    // distance to group i.  Replaces the O(ntrans) scan in merge() with an
    // O(degree) traversal — decisive for the 500k-transcript super-cluster.
    std::vector<std::vector<int>> adj_;
    std::vector<bool>             alive_;  // false once a group is merged away

    int id_ = 0;
    std::vector<int> sample_groups;

    // --- Private clustering methods ---
    unsigned char find_next_pair(int &max_i, int &max_j);
    void          merge(int i, int j);
    void          initialise_matrix();
    std::vector<int> get_counts(int s);
    void          output_clusters(const std::string &threshold,
                                  const std::string &method_tag = "");

public:
    // --- Accessors ---
    void        add_tran(Transcript *t) { cluster_.push_back(t); }
    void        add_read(Read *r)       { read_.push_back(r); }
    [[nodiscard]] Transcript *get_tran(int i) const { return cluster_[i]; }
    [[nodiscard]] Read       *get_read(int i) const { return read_[i]; }
    [[nodiscard]] int         n_trans() const       { return static_cast<int>(cluster_.size()); }
    [[nodiscard]] int         n_reads() const       { return static_cast<int>(read_.size()); }
    [[nodiscard]] int         get_id() const        { return id_; }
    void        set_id(int id)   { id_ = id; }

    // Index into MakeClusters::clusterList for O(1) swap-and-pop removal
    int clusterList_index = -1;

    void set_sample_groups(const std::vector<int> &sg) { sample_groups = sg; }

    // Distance between two transcript groups (shared-read proportion
    // with LRT expression test).  Public so LeidenCluster can reuse it.
    // When lrt_softness > 0, applies a sigmoid decay around D_cut instead
    // of a hard binary gate (Leiden-only; hierarchical uses default 0).
    float get_dist(int i, int j, float lrt_softness = 0.0f);

    // Build read_groups / read_group_sizes from reads.  Called by
    // initialise_matrix() and cluster_leiden() before graph construction.
    void setup_read_groups();

    // Accessors for Leiden graph construction
    [[nodiscard]] const group &get_groups() const         { return groups; }
    [[nodiscard]] const read_group &get_read_groups() const { return read_groups; }
    [[nodiscard]] const group &get_read_group_sizes() const { return read_group_sizes; }

    // Main entry point: hierarchical clustering with expression testing.
    void cluster(std::map<float, std::string> &thresholds,
                 const std::string &method_tag = "");

    // --- Static configuration ---
    static inline ClusterAlgorithm algorithm = ClusterAlgorithm::Hierarchical;
    static inline float D_cut = 0;
    static inline float lrt_softness = 0.0f;  // sigmoid steepness (0 = hard cutoff)
    static inline int   knn = 0;              // -1 = disabled, 0 = auto (default), >0 = fixed k
    static inline std::string file_prefix;
    static constexpr std::string_view file_counts              = "counts";
    static constexpr std::string_view file_clusters            = "clusters";
    static constexpr std::string_view file_ext                 = ".txt";
    static constexpr std::string_view cluster_id_prefix        = "Cluster-";
    static constexpr std::string_view cluster_id_joiner        = ".";
    static constexpr std::string_view cluster_id_prefix_no_reads = "NoReadsCluster-";

    void print_alignments() const;
};
