// Leiden community detection wrapper using the igraph C library.
// Provides an alternative to the O(n^2) hierarchical clustering
// in Cluster.cc.  Uses the Constant Potts Model (CPM) quality
// function to avoid the resolution-limit problem.
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
// Cluster::cluster() when --algorithm leiden is selected.
void cluster_leiden(Cluster *c,
                    std::map<float, std::string> &thresholds);

#endif // HAVE_IGRAPH
