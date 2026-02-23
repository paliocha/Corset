// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.
//
// Last modified 21 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#include <Transcript.h>
#include <Read.h>

void Transcript::remove() {
    for (auto *r : reads_)
        r->remove(this);
}

void Transcript::add_read(Read *read) {
    reads_.push_back(read);
}

bool Transcript::reached_min_counts() const {
    int counts = total_direct_counts();
    if (counts >= min_counts)
        return true;
    for (auto *r : reads_) {
        counts += r->get_weight();
        if (counts >= min_counts)
            return true;
    }
    return false;
}

