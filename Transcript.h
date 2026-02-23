// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// Stores transcript names and temporary read-back-references.
// Allows transcripts to be identified by pointer rather than
// string, which is faster and uses less memory.  Also holds
// global clustering parameters as static members (samples,
// groups, min_counts) set once from CLI arguments.
//
// Original author: Nadia Davidson
// Last modified 21 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#pragma once

#include <string>
#include <vector>
#include <StringSet.h>

class Read;

class Transcript {
    std::string name_;
    int pos_ = 0;                  // position index set by Cluster::initialise_matrix()
    std::vector<Read *> reads_;    // back-references for min-count filtering
    std::vector<int> direct_counts_;  // [sample] — lazy-init, for redistributed ECs

public:
    Transcript()                          : name_("") {}
    explicit Transcript(const std::string &name) : name_(name) {}

    [[nodiscard]] const std::string &get_name() const { return name_; }
    void  pos(int position)             { pos_ = position; }
    [[nodiscard]] int pos() const       { return pos_; }

    void add_read(Read *read);
    [[nodiscard]] bool reached_min_counts() const;
    void remove();   // remove this transcript from its reads' alignment lists

    std::vector<Read *> *get_reads() { return &reads_; }

    // Direct counts — per-sample counters for redistributed ECs (weight < -l).
    // Lazy-initialized: only transcripts that receive redistributed reads allocate.
    void add_direct_count(int sample, int weight) {
        if (direct_counts_.empty())
            direct_counts_.resize(samples, 0);
        direct_counts_[sample] += weight;
    }
    [[nodiscard]] int get_direct_count(int sample) const {
        return direct_counts_.empty() ? 0 : direct_counts_[sample];
    }
    [[nodiscard]] int total_direct_counts() const {
        int sum = 0;
        for (int c : direct_counts_) sum += c;
        return sum;
    }
    [[nodiscard]] bool has_direct_counts() const { return !direct_counts_.empty(); }

    // Global configuration (set once from CLI args)
    static inline int samples          = 0;
    static inline int groups           = 0;
    static inline int min_counts       = 10;
    static inline int min_reads_for_link = 1;
    static inline int max_alignments   = -1;
};

using TranscriptList = StringSet<Transcript>;
