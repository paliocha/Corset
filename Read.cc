// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.
//
// Last modified 21 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#include <Read.h>
#include <unordered_map>

using std::unordered_map;
using std::vector;
using std::string;
using std::ofstream;
using std::endl;
using std::cout;

// BAM path: look up read + transcript by name
void ReadList::add_alignment(const string &read, const string &trans, int sample) {
    Read *r       = reads_map->insert(read);
    Transcript *t = transcript_list->insert(trans);
    if (!r->has(t)) {
        r->add_alignment(t);
        r->set_sample(sample);
    }
}

// Corset / string eq_class path: transcript names resolved via StringSet
void ReadList::add_alignment(vector<string> &trans_names, int sample, int weight) {
    Read *r = new Read();
    reads_vector.push_back(r);
    r->set_sample(sample);
    r->set_weight(weight);
    for (auto &name : trans_names) {
        Transcript *t = transcript_list->insert(name);
        r->add_alignment(t);
    }
}

// Salmon eq_class fast path: Transcript pointers already resolved
void ReadList::add_alignment_ptrs(vector<Transcript *> &trans_ptrs, int sample, int weight) {
    Read *r = new Read();
    reads_vector.push_back(r);
    r->set_sample(sample);
    r->set_weight(weight);
    for (Transcript *t : trans_ptrs)
        r->add_alignment(t);
}

// Merge reads with identical alignment sets into weighted compact reads.
// Uses hash-bucket grouping to avoid O(n^2) pairwise comparison.
void ReadList::compactify_reads(TranscriptList *trans, const string & /*outputReadsName*/) {
    // Sort alignments and compute hashes for all reads
    for (auto &[name, read] : *reads_map) {
        read->sort_alignments();
        read->set_trans_hash();
    }

    // For each transcript, group its reads by hash and merge duplicates
    for (auto &[name, transcript] : *trans) {
        vector<Read *> *reads = transcript->get_reads();
        int reads_size = static_cast<int>(reads->size());
        if (reads_size < 2) continue;

        // Build hash → indices map
        unordered_map<uintptr_t, vector<int>> hash_buckets;
        hash_buckets.reserve(reads_size);
        for (int i = 0; i < reads_size; i++) {
            if ((*reads)[i]->get_weight() != 0)
                hash_buckets[(*reads)[i]->get_trans_hash()].push_back(i);
        }

        // Only compare within each hash bucket
        for (auto &[hash, indices] : hash_buckets) {
            for (int bi = 0; bi < static_cast<int>(indices.size()) - 1; bi++) {
                Read *ri = (*reads)[indices[bi]];
                if (ri->get_weight() == 0) continue;
                for (int bj = bi + 1; bj < static_cast<int>(indices.size()); bj++) {
                    Read *rj = (*reads)[indices[bj]];
                    if (rj->get_weight() != 0 && ri->has_same_alignments(rj)) {
                        ri->set_weight(ri->get_weight() + rj->get_weight());
                        rj->set_weight(0);
                    }
                }
            }
        }

        // Compact: remove zero-weight reads from this transcript's list
        int write_pos = 0;
        for (int k = 0; k < static_cast<int>(reads->size()); k++) {
            if ((*reads)[k]->get_weight() != 0)
                (*reads)[write_pos++] = (*reads)[k];
        }
        reads->resize(write_pos);
    }

    // Move surviving reads into the flat vector; delete zero-weight reads
    for (auto &[name, r] : *reads_map) {
        if (r->get_weight() != 0)
            reads_vector.push_back(r);
        else
            delete r;
    }
    reads_map->clear();
    delete reads_map;
    reads_map = nullptr;
}

void ReadList::write(const string &outputReadsName) {
    ofstream readFile(outputReadsName);
    for (Read *r : reads_vector) {
        readFile << r->get_weight();
        for (auto trans = r->align_begin(); trans != r->align_end(); ++trans)
            readFile << "\t" << (*trans)->get_name();
        readFile << "\n";
    }
    readFile.close();
    cout << "Done writing " << outputReadsName << endl;
}
