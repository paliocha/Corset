# Changelog — corset-omp

## 1.10.0 (2026-02-20)

Fork of [Oshlack/Corset](https://github.com/Oshlack/Corset) v1.09 with
OpenMP parallelism, modern htslib, and a full C++23 port.

### OpenMP parallelism
- Hierarchical clustering dispatched in parallel across super-clusters
  (`#pragma omp parallel for schedule(dynamic)`).
- New CLI flag `-t <int>` to set thread count (default: `OMP_NUM_THREADS`).
- Thread-safe file I/O via `#pragma omp critical` with persistent handles.
- Thread-safe RNG (`rand_r`) in `get_counts()`.

### htslib 1.23 port
- Replaced legacy samtools 0.x inline API with htslib 1.23 (`sam_open`,
  `sam_hdr_read`, `sam_read1`, `bam_get_qname`, `sam_hdr_tid2name`).
- Container builds htslib 1.23 from source with libcurl support.

### Performance optimisations
- **Adjacency-list merge** (`adj_`/`alive_`): Replace O(ntrans) linear scan
  in `merge()` with O(degree) sorted-merge of per-group neighbor lists.
  For the 543k-transcript super-cluster, each merge visits ~50–500 live
  neighbors instead of scanning all 543k indices.
- **Parallel distance recomputation**: After collecting merged neighbors,
  `get_dist()` calls are parallelized with `#pragma omp parallel for`
  (nested, threshold ncand > 128).  `get_dist()` is read-only on shared
  data, safe for concurrent execution.
- **Flat open-addressing hash map** (`FlatDistMap.h`): Custom linear-probing
  hash map with Wang hash replacing `std::unordered_map` for the distance
  matrix.  Flat memory layout, no per-node allocation, ~2–3× better cache
  performance for point lookups.
- **Early-out per sample** in `get_dist()`: Skip sorted-list intersection
  entirely when either transcript group has zero reads in a sample.
  ~30–50% of samples are zero for any given pair in typical RNA-Seq data.
- **SSE2 block intersection**: Process sorted list A in 4-element blocks,
  compare all qualifying B elements via `_mm_cmpeq_epi32` (4 comparisons
  per instruction).  Falls through to scalar merge for remainder.
- **Lazy-deletion max-heap** for `find_next_pair()` — avoids O(n²) scan.
- **Hash-bucket read deduplication** in `compactify_reads()` — avoids O(n²)
  pairwise comparison.
- **Weighted union-find** for super-cluster construction — O(n log n)
  transcript reassignment.
- **Pre-resolved transcript pointers** in salmon eq-class reader — O(1)
  integer indexing instead of string lookups per equivalence class.
- **Incremental read-group size tracking** during merge — avoids re-scanning
  all reads.

### C++23 standard port (from C++0x/C++11)
Build system:
- `configure.in`: `-std=c++0x -DUNORDEREDMAP` → `-std=c++23`.
- `Makefile.in`: `corset_fasta_ID_changer` now built with `@CXXFLAGS@`.
- Container base: Ubuntu 22.04 → 24.04 (GCC 13.2, full C++23 support).

Headers & language features:
- `#pragma once` on all headers (replacing include guards).
- `[[nodiscard]]` on all pure getters and query methods.
- `static inline` member initialisers (eliminating .cc boilerplate).
- `static constexpr string_view` for all string constants — zero heap allocations.
- `static constexpr operator()` on `WangHash64` (C++23) — no `this` pointer.
- `constexpr` on `HeapEntry::operator<`.
- `enum class ReadStatus` replacing `#define` constants.

Containers & algorithms:
- `unordered_map::contains()` replacing `find() == end()` (C++20).
- `std::ranges::sort` / `std::ranges::unique` (C++20) in 5 call sites.
- `std::ranges::contains()` replacing `std::find() != end()` (C++23).
- `std::erase()` replacing find+erase pattern (C++20).
- `try_emplace()` in `StringSet::insert` — single hash lookup (C++17).
- `DistanceMatrix::get()` made `const` and find-based — prevents spurious
  default-insertions into the sparse map.

String operations:
- `string::contains()` replacing `find() != npos` (C++23).
- `string::starts_with()` replacing `find('>') == 0` (C++20).
- `std::stoi` replacing C `atoi` for safer integer parsing.

Modern idioms:
- Structured bindings in all map/pair iterations.
- Range-for loops replacing iterator boilerplate throughout.
- VLAs replaced with `std::vector` (standards compliance).
- Explicit `#include <cstdint>` for `uint64_t` / `uintptr_t` — required
  under strict C++23 (GCC 13 no longer provides these transitively).

### Code hygiene (carried from pre-port cleanup)
- Removed all `using namespace std` directives.
- Consistent formatting, dead code removal, const-correctness.
- Modern `using` aliases replacing `typedef`.
