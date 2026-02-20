// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// A Cluster (super-cluster) holds transcripts that share at least one
// read.  Contains the hierarchical clustering algorithm and the
// log-likelihood ratio expression test.
//
// Original author: Nadia Davidson

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
#include <unordered_map>
#include <cstdint>
#include <climits>
#include <Read.h>

using group       = std::vector<std::vector<int>>;
using read_group  = std::vector<std::vector<std::vector<int>>>;

// Thomas Wang 64-bit integer hash — excellent distribution for
// structured keys like ntrans*i+j used in the sparse distance matrix.
struct WangHash64 {
    [[nodiscard]] static constexpr size_t operator()(uint64_t key) noexcept {
        key = (~key) + (key << 21);
        key ^= (key >> 24);
        key = (key + (key << 3)) + (key << 8);
        key ^= (key >> 14);
        key = (key + (key << 2)) + (key << 4);
        key ^= (key >> 28);
        key += (key << 31);
        return static_cast<size_t>(key);
    }
};

using dist_map      = std::unordered_map<uint64_t, unsigned char, WangHash64>;
using dist_iterator = dist_map::iterator;

// Lazy-deletion max-heap entry for find_next_pair().
struct HeapEntry {
    unsigned char value;
    int i, j;   // invariant: i > j
    [[nodiscard]] constexpr bool operator<(const HeapEntry &o) const { return value < o.value; }
};

// Sparse distance matrix keyed by (i*ntrans + j).
class DistanceMatrix {
    dist_map dist_;
    uint64_t ntrans_ = 0;

public:
    void set_size(int n)    { ntrans_ = static_cast<uint64_t>(n); }
    void reserve(size_t n)  { dist_.reserve(n); }

    [[nodiscard]] unsigned char get(int i, int j) const {
        auto it = dist_.find(ntrans_ * i + j);
        return it != dist_.end() ? it->second : 0;
    }
    void set(int i, int j, int value) {
        if (value != 0) dist_[ntrans_ * i + j] = static_cast<unsigned char>(value);
        else            remove(i, j);
    }
    [[nodiscard]] bool no_link(int i, int j) const { return !dist_.contains(ntrans_ * i + j); }
    void remove(int i, int j)                      { dist_.erase(ntrans_ * i + j); }

    dist_iterator begin() { return dist_.begin(); }
    dist_iterator end()   { return dist_.end(); }
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

    int id_ = 0;
    std::vector<int> sample_groups;

    // --- Private clustering methods ---
    float         get_dist(int i, int j);
    unsigned char find_next_pair(int &max_i, int &max_j);
    void          merge(int i, int j);
    void          initialise_matrix();
    std::vector<int> get_counts(int s);
    void          output_clusters(const std::string &threshold);

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

    // Main entry point: hierarchical clustering with expression testing.
    void cluster(std::map<float, std::string> &thresholds);

    // --- Static configuration ---
    static inline float D_cut = 0;
    static inline std::string file_prefix;
    static constexpr std::string_view file_counts              = "counts";
    static constexpr std::string_view file_clusters            = "clusters";
    static constexpr std::string_view file_ext                 = ".txt";
    static constexpr std::string_view cluster_id_prefix        = "Cluster-";
    static constexpr std::string_view cluster_id_joiner        = ".";
    static constexpr std::string_view cluster_id_prefix_no_reads = "NoReadsCluster-";

    void print_alignments() const;
};
