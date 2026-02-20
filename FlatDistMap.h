// Cache-friendly open-addressing hash map with linear probing.
// Keys: uint64_t, Values: unsigned char.
//
// ~2-3× faster than std::unordered_map for point operations: no pointer
// chasing, no per-node heap allocation, better prefetch behavior.
//
// Uses sentinel keys (EMPTY/TOMBSTONE) instead of a separate metadata
// array.  Safe key range: [0, UINT64_MAX-2].  The distance-matrix keys
// (ntrans*i + j, ntrans ≤ 2M, i,j < ntrans) are well within range.

#pragma once

#include <cstdint>
#include <vector>

class FlatDistMap {
    static constexpr uint64_t EMPTY     = UINT64_MAX;
    static constexpr uint64_t TOMBSTONE = UINT64_MAX - 1;
    static constexpr double   MAX_LOAD  = 0.7;

    struct Slot {
        uint64_t      key   = EMPTY;
        unsigned char value = 0;
    };

    std::vector<Slot> slots_;
    size_t size_     = 0;
    size_t capacity_ = 0;   // always power of 2
    size_t mask_     = 0;   // capacity_ - 1

    // WangHash64 inlined — avoids indirection through a functor.
    [[nodiscard]] static constexpr size_t hash(uint64_t k) noexcept {
        k = (~k) + (k << 21);
        k ^= (k >> 24);
        k = (k + (k << 3)) + (k << 8);
        k ^= (k >> 14);
        k = (k + (k << 2)) + (k << 4);
        k ^= (k >> 28);
        k += (k << 31);
        return static_cast<size_t>(k);
    }

    void insert_no_grow(uint64_t key, unsigned char value) {
        size_t idx = hash(key) & mask_;
        while (true) {
            auto &s = slots_[idx];
            if (s.key == EMPTY || s.key == TOMBSTONE) {
                s = {key, value};
                ++size_;
                return;
            }
            if (s.key == key) { s.value = value; return; }
            idx = (idx + 1) & mask_;
        }
    }

    void grow_and_rehash() {
        size_t new_cap = capacity_ ? capacity_ * 2 : 1024;
        std::vector<Slot> old = std::move(slots_);
        slots_.assign(new_cap, Slot{});
        capacity_ = new_cap;
        mask_     = new_cap - 1;
        size_     = 0;
        for (const auto &s : old)
            if (s.key != EMPTY && s.key != TOMBSTONE)
                insert_no_grow(s.key, s.value);
    }

public:
    FlatDistMap() = default;

    void reserve(size_t n) {
        size_t needed = static_cast<size_t>(n / MAX_LOAD) + 1;
        size_t cap = 1024;
        while (cap < needed) cap *= 2;
        if (cap > capacity_) {
            std::vector<Slot> old = std::move(slots_);
            slots_.assign(cap, Slot{});
            capacity_ = cap;
            mask_     = cap - 1;
            size_     = 0;
            for (const auto &s : old)
                if (s.key != EMPTY && s.key != TOMBSTONE)
                    insert_no_grow(s.key, s.value);
        }
    }

    void set(uint64_t key, unsigned char value) {
        if (size_ >= static_cast<size_t>(capacity_ * MAX_LOAD))
            grow_and_rehash();
        insert_no_grow(key, value);
    }

    [[nodiscard]] unsigned char get(uint64_t key) const noexcept {
        if (capacity_ == 0) return 0;
        size_t idx = hash(key) & mask_;
        while (true) {
            const auto &s = slots_[idx];
            if (s.key == key)   return s.value;
            if (s.key == EMPTY) return 0;
            idx = (idx + 1) & mask_;
        }
    }

    [[nodiscard]] bool contains(uint64_t key) const noexcept {
        if (capacity_ == 0) return false;
        size_t idx = hash(key) & mask_;
        while (true) {
            if (slots_[idx].key == key)   return true;
            if (slots_[idx].key == EMPTY) return false;
            idx = (idx + 1) & mask_;
        }
    }

    void erase(uint64_t key) noexcept {
        if (capacity_ == 0) return;
        size_t idx = hash(key) & mask_;
        while (true) {
            auto &s = slots_[idx];
            if (s.key == key) { s.key = TOMBSTONE; --size_; return; }
            if (s.key == EMPTY) return;
            idx = (idx + 1) & mask_;
        }
    }

    [[nodiscard]] size_t size() const noexcept { return size_; }
};
