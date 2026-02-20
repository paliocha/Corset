// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// Read and ReadList classes.
// A Read records the transcript alignments for one (possibly weighted)
// equivalence class, tagged with a sample index.
//
// Original author: Nadia Davidson

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <Transcript.h>

// A single read (or weighted equivalence class).
class Read {
    std::vector<Transcript *> alignments_;
    unsigned char sample_ = 0;   // max 255 samples (see audit #15)
    int weight_           = 1;
    uintptr_t trans_hash_ = 0;

public:
    Read() = default;
    explicit Read(const std::string & /*name*/) {}  // dummy ctor for StringSet

    void set_sample(int s)                     { sample_ = static_cast<unsigned char>(s); }
    [[nodiscard]] int  get_sample() const      { return sample_; }
    [[nodiscard]] int  alignments() const      { return static_cast<int>(alignments_.size()); }
    void set_weight(int w)                     { weight_ = w; }
    [[nodiscard]] int  get_weight() const      { return weight_; }

    void add_alignment(Transcript *t) {
        alignments_.push_back(t);
        t->add_read(this);
    }

    using trans_iter = std::vector<Transcript *>::iterator;
    trans_iter align_begin() { return alignments_.begin(); }
    trans_iter align_end()   { return alignments_.end(); }

    void sort_alignments() { std::ranges::sort(alignments_); }

    [[nodiscard]] bool has(Transcript *t) const {
        return std::ranges::contains(alignments_, t);
    }

    void remove(Transcript *t) { std::erase(alignments_, t); }

    // Check if two reads align to exactly the same transcripts.
    // Requires sort_alignments() + set_trans_hash() to have been called first.
    [[nodiscard]] bool has_same_alignments(const Read *r) const {
        return trans_hash_ == r->trans_hash_ && sample_ == r->sample_;
    }

    [[nodiscard]] uintptr_t get_trans_hash() const { return trans_hash_; }

    // Compute a combined hash over aligned transcript pointers (boost::hash_combine style).
    void set_trans_hash() {
        trans_hash_ = 0;
        for (auto *t : alignments_)
            trans_hash_ ^= std::hash<uintptr_t>()(reinterpret_cast<uintptr_t>(t))
                           + 0x9e3779b9 + (trans_hash_ << 6) + (trans_hash_ >> 2);
    }

    void print_alignments() const {
        for (auto *t : alignments_)
            std::cout << t->get_name() << " ";
        std::cout << std::endl;
    }
};


// A container for a set of Reads belonging to one sample.
class ReadList {
    TranscriptList *transcript_list;
    StringSet<Read> *reads_map;          // destroyed after compactify_reads()
    std::vector<Read *> reads_vector;

public:
    explicit ReadList(TranscriptList *transcripts)
        : transcript_list(transcripts), reads_map(new StringSet<Read>) {}

    // BAM input path — looks up read name + transcript name via StringSet
    void add_alignment(const std::string &read, const std::string &trans, int sample);

    // Corset / string-based eq_class input
    void add_alignment(std::vector<std::string> &trans_names, int sample, int weight);

    // Salmon eq_class fast path — pre-resolved Transcript pointers
    void add_alignment_ptrs(std::vector<Transcript *> &trans_ptrs, int sample, int weight);

    // Merge duplicate reads (same alignment set) into weighted compact reads.
    void compactify_reads(TranscriptList *trans, const std::string &outputReadsName = "");

    std::vector<Read *>::iterator begin() { return reads_vector.begin(); }
    std::vector<Read *>::iterator end()   { return reads_vector.end(); }

    void write(const std::string &outputReadsName);
};
