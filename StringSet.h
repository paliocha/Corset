// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// A template container storing a pointer to an object
// alongside its name (ID string).  Used for Read and
// Transcript lookup by name during I/O; destroyed once
// super-clusters are built to conserve memory.
//
// Original author: Nadia Davidson
// Last modified 21 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#pragma once

#include <string>
#include <cstdlib>
#include <unordered_map>

template <class T>
class StringSet {
    using map_type = std::unordered_map<std::string, T *>;

    map_type set_map;

public:
    using iterator = typename map_type::iterator;

    // Look up an object by name.  Returns nullptr if not found.
    [[nodiscard]] T *find(const std::string &name) {
        auto it = set_map.find(name);
        return (it != set_map.end()) ? it->second : nullptr;
    }

    // Insert a new object keyed by name.  If the name already exists,
    // return the existing object pointer.
    // Uses try_emplace for a single hash lookup (find + insert combined).
    T *insert(const std::string &name) {
        auto [it, inserted] = set_map.try_emplace(name, nullptr);
        if (!inserted)
            return it->second;
        it->second = new T(name);
        return it->second;
    }

    iterator begin() { return set_map.begin(); }
    iterator end()   { return set_map.end(); }
    void clear()     { set_map.clear(); }
    [[nodiscard]] int size() const { return static_cast<int>(set_map.size()); }
};
