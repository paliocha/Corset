# LeidenCluster.cc — Post-audit improvements

Issues from code audit, ordered for implementation (dependencies first,
then descending priority).

---

## 1. Negative Q_ref comparison inversion (correctness)

**Issue #5.** CPM quality can be negative at high resolution. The ratio
`knn_result.quality / Q_ref` inverts when Q_ref < 0: a worse (more
negative) quality scores > 1.0 and passes the threshold.

**File:** [LeidenCluster.cc:326-327](LeidenCluster.cc#L326-L327)

**Fix:** Replace ratio with direct threshold comparison:
```cpp
// Before:
double quality_ratio = (std::abs(Q_ref) > 1e-15)
    ? knn_result.quality / Q_ref : 1.0;
bool ok = (quality_ratio >= KNN_QUALITY_THRESHOLD) && ...

// After:
bool quality_ok = (std::abs(Q_ref) > 1e-15)
    ? knn_result.quality >= KNN_QUALITY_THRESHOLD * Q_ref
    : true;
bool ok = quality_ok && (nmi >= KNN_NMI_THRESHOLD);
```

Works for positive Q_ref (0.8 * Q scales down), negative Q_ref (0.8 *
negative is less negative = easier threshold), and zero Q_ref (skip).

---

## 2. Zero-edge short circuit in kNN auto path (correctness/clarity)

**Issue #7.** When all pairs have `get_dist() == 0`, the reference Leiden
run gets an edgeless graph. The binary search runs ~6 probes on empty
edge lists, all returning empty. Correct but wasteful, and the log line
prints confusing "0 -> 0 edges".

**File:** [LeidenCluster.cc:528-556](LeidenCluster.cc#L528-L556) (kNN
auto block, after `collect_edges_knn`)

**Fix:** Early return after edge collection:
```cpp
edges = collect_edges_knn(c, ntrans, k_max);
if (edges.size() == 0) {
    // All pairs LRT-zeroed — skip kNN search, fall through to
    // per-threshold Leiden which will produce singletons.
    goto build_graph;
}
```
Or restructure with a flag / nested block to avoid goto.

---

## 3. `node_adj` memory pressure in Phase 1 (memory)

**Issue #1.** `vector<vector<int>> node_adj(ntrans)` at
[LeidenCluster.cc:117](LeidenCluster.cc#L117) stores the full adjacency
list for every node. For a dense 50k-transcript SC with avg degree 5000,
that's ~1 GB concurrent with `pair_counts`.

**Fix options (pick one):**

**(a) Free each node's adjacency after Phase 1b processes it.**
Phase 1b iterates nodes 0..ntrans-1 sequentially. After selecting top-k
neighbors for node `n`, do `node_adj[n] = {}; node_adj[n].shrink_to_fit();`.
Requires Phase 1 to finish completely before Phase 1b starts (already the
case). Reduces peak from O(total_adjacency) to O(max_single_node_degree).

**(b) Add iteration to FlatDistMap, eliminate `node_adj` entirely.**
Add a `for_each(callback)` method that scans the slot array. Phase 1
only populates `pair_counts`. Phase 1b iterates `pair_counts`, decodes
keys back to (hi, lo) via `hi = key / ntrans, lo = key % ntrans`, and
builds `scored` per-node on the fly. Trades one pass over the hash table
for eliminating O(edges) adjacency storage. More invasive but cleaner.

Option (a) is simpler and sufficient.

---

## 4. Per-node `scored` vector allocation in Phase 1b (performance)

**Issue #3.** [LeidenCluster.cc:149-161](LeidenCluster.cc#L149-L161)
allocates a `vector<pair<unsigned char, int>>` for each node with > k_max
neighbors, totaling up to 50k heap allocations.

**Fix:** Hoist a reusable buffer outside the node loop:
```cpp
vector<std::pair<unsigned char, int>> scored;
for (int n = 0; n < ntrans; n++) {
    auto &adj = node_adj[n];
    if (static_cast<int>(adj.size()) <= k_max) { ... continue; }

    scored.clear();
    scored.reserve(adj.size());   // no-op after first large node
    for (int nb : adj) { ... }
    std::nth_element(...);
    ...
}
```

Combine with issue #3 fix: free `node_adj[n]` after processing.

---

## 5. `selected` reserve underestimate (performance)

**Issue #2.** [LeidenCluster.cc:139](LeidenCluster.cc#L139) reserves
`ntrans * min(k_max, 64)` slots, but actual selected pairs can reach
`ntrans * k_max / 2`. FlatDistMap rehashes mid-insertion when k_max > 64.

**Fix:**
```cpp
// Before:
selected.reserve(static_cast<size_t>(ntrans) * std::min(k_max, 64));

// After:
selected.reserve(static_cast<size_t>(ntrans) * k_max / 2);
```

The `/ 2` accounts for symmetric dedup (each undirected edge is stored
once under canonical key hi*ntrans+lo, but nominated by up to 2 nodes).

---

## 6. Phase 2 iterates full adjacency to check `selected` (performance)

**Issue #4.** [LeidenCluster.cc:172-184](LeidenCluster.cc#L172-L184)
loops over every neighbor in `node_adj[n]` and does a hash lookup in
`selected` for each. For a node with 5000 neighbors and k_max=224,
that's ~5000 lookups to emit ~224 edges.

**Fix options:**

**(a) Iterate `selected` instead of `node_adj`.**
Requires FlatDistMap iteration (see issue #3 option b). Scan `selected`,
decode each key to (hi, lo), call `get_dist(hi, lo)`, emit if > 0. Cost:
O(selected.capacity) scan vs O(total_adjacency) with hash lookups. Better
when selection rate is low.

**(b) Keep current approach but free `node_adj[n]` in Phase 1b.**
If issue #3 option (a) is chosen, `node_adj` is already freed by the time
Phase 2 runs. Phase 2 would need a different source of candidate pairs.
This forces option (a) here too — iterate `selected`.

**(c) Accept current cost.** The hash lookups are O(1) each, and Phase 2's
bottleneck is `get_dist()` (SSE2 sorted-list intersection), not the lookup.
The adjacency iteration adds ~5% overhead vs the distance computation.
Pragmatically fine.

Recommendation: option (c) for now. Revisit if profiling shows otherwise.

---

## 7. `collect_edges` reserve heuristic (performance, minor)

**Issue #10.** [LeidenCluster.cc:73](LeidenCluster.cc#L73) reserves
`ntrans * 4` for the seen-set. Underestimates for dense SCs, causing
rehashes.

**Fix:** Use `ntrans * k_avg_estimate` where `k_avg_estimate` is derived
from the read count: `min(c->n_reads(), ntrans * ntrans / 2)` capped at
a reasonable maximum. Or simply accept the current heuristic — FlatDistMap
rehash is O(n) amortized and only happens on the non-kNN path (kNN path
uses `collect_edges_knn` which has its own reserve).

No action needed unless profiling shows rehash as a bottleneck.

---

## 8. `apply_knn` copies EdgeList per binary search probe (memory)

**Issue #6.** [LeidenCluster.cc:312](LeidenCluster.cc#L312) calls
`apply_knn(edges, ntrans, mid)` which returns a new EdgeList each
iteration. For ~6 probes on a 10k SC with 500k edges, that's 6 × 12 MB
= 72 MB of transient allocation.

**Fix:** Not worth optimizing — the copies are short-lived (destroyed
at end of each loop iteration) and small relative to the graph
construction inside `build_igraph`. The binary search only runs ~6 times.

No action needed.

---

## 9. `rand_r` determinism note (correctness, non-issue)

**Issue #8.** [LeidenCluster.cc:397,424](LeidenCluster.cc#L397) uses
`rand_r` with seed `id * 17 + sample * 31 + 42`. Thread-safe, seed is
thread-independent. Reproducible across runs regardless of thread count.

No action needed.

---

## 10. igraph global RNG data race (correctness, upstream limitation)

**Issue #9.** [LeidenCluster.cc:515-516](LeidenCluster.cc#L515-L516)
calls `igraph_rng_seed(igraph_rng_default(), ...)` on the global default
RNG from inside an OMP parallel region. `igraph_rng_default()` is not
thread-safe — concurrent seed/draw calls are a data race.

**Impact:** Low. Leiden CPM on weighted graphs converges to the same
partition regardless of initial assignment in practice. The RNG affects
only the random node-visit order in refinement steps.

**Fix options:**

**(a) Remove the seed call entirely.** Let igraph use its default seed.
Leiden results are practically deterministic for CPM anyway.

**(b) Seed once before the OMP parallel region.** Loses per-SC seeding
but eliminates the race. Sufficient if exact per-SC reproducibility is
not required.

**(c) Per-thread igraph RNG.** `igraph_community_leiden_simple()` does
not accept a custom RNG — would require switching to the lower-level
`igraph_community_leiden()` API (if available in igraph 1.0.1). Most
correct, most invasive.

Recommendation: option (a) — delete the two lines. The seeding provides
negligible value and introduces a race.

---

## Implementation order

| Step | Issue | Type | Effort |
|------|-------|------|--------|
| 1 | #1 — Negative Q_ref | correctness | 5 min |
| 2 | #2 — Zero-edge short circuit | correctness | 5 min |
| 3 | #10 — Remove igraph RNG seed | data race | 1 min |
| 4 | #3 — Free `node_adj` incrementally | memory | 10 min |
| 5 | #4 — Hoist `scored` buffer | performance | 5 min |
| 6 | #5 — Fix `selected` reserve | performance | 1 min |
| 7-10 | #6-#9 | no action | — |

Total: ~30 minutes for issues 1-6. Issues 7-10 require no changes.
