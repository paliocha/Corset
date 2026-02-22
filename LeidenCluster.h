// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
//
// Features:
//   - Soft LRT weighting (--lrt-softness): sigmoid decay around D_cut
//   - Adaptive kNN sparsification (--knn): connectivity-constrained
//     minimum-k selection per super-cluster
//
// Author: Martin Paliocha, martin.paliocha@nmbu.no
// Created 22 February 2026

#pragma once

#ifdef HAVE_IGRAPH

#include <map>
#include <string>

class Cluster;

// Run Leiden CPM on one super-cluster, write output files.
// Called from MakeClusters::processSuperClusters() in place of
// Cluster::cluster() when --algorithm leiden or both is selected.
// method_tag: "" for standalone, "l-" for comparison mode.
void cluster_leiden(Cluster *c,
                    std::map<float, std::string> &thresholds,
                    const std::string &method_tag = "");

#endif // HAVE_IGRAPH
