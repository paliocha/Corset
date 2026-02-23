// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
//
// Features:
//   - Soft LRT weighting (--lrt-softness): sigmoid decay around D_cut
//   - Adaptive kNN sparsification (--knn): quality-preserving
//     minimum-k selection per super-cluster (CPM quality + NMI)
//   - Two-phase edge construction: cheap proxy ranking + get_dist
//     only on kNN candidates (avoids O(E) distance computation)
//
// Author: Martin Paliocha, martin.paliocha@nmbu.no
// Created 22 February 2026

#pragma once

#ifdef HAVE_IGRAPH

#include <map>
#include <string>
#include <vector>

class Cluster;

// Run Leiden CPM on one super-cluster, write output files.
// Called from MakeClusters::processSuperClusters() in place of
// Cluster::cluster() when --algorithm leiden or both is selected.
// method_tag: "" for standalone, "l-" for comparison mode.
// cluster_sizes: if non-null, filled with per-community transcript counts
//                from the last threshold (for summary statistics).
// Returns the number of edges in the (possibly kNN-filtered) graph.
int cluster_leiden(Cluster *c,
                   std::map<float, std::string> &thresholds,
                   const std::string &method_tag = "",
                   std::vector<int> *cluster_sizes = nullptr);

#endif // HAVE_IGRAPH
