// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// Builds super-clusters (stage 1) via weighted union-find,
// then dispatches OpenMP-parallel clustering (stage 2) across
// independent super-clusters.  Supports hierarchical (default)
// and Leiden (--algorithm leiden) back-ends.
//
// Original author: Nadia Davidson
// Last modified 23 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <Cluster.h>
#include <Read.h>
#include <Transcript.h>

// Per-super-cluster timing for the summary report.
struct SCStats {
    int id;
    int ntrans;
    int nedges;
    int nclusters;
    double duration;
};

class MakeClusters {
    std::vector<Cluster *> clusterList;
    Cluster *current_cluster = nullptr;

    std::unordered_map<Transcript *, Cluster *> transMap;

    // Find-or-create the cluster for a transcript
    std::pair<Transcript *const, Cluster *> *getMapElement(Transcript *trans);

    void setCurrentCluster(Transcript *trans) {
        current_cluster = getMapElement(trans)->second;
    }
    void checkAgainstCurrentCluster(Transcript *trans);

    void makeSuperClusters(std::vector<ReadList *> &readLists,
                           TranscriptList *tList);
    void processSuperClusters(std::map<float, std::string> &thresholds,
                              std::vector<int> &groups);

    // Summary data (populated during stages, printed at the end)
    double time_reading_        = 0;
    double time_superclusters_  = 0;
    double time_clustering_     = 0;
    int    total_transcripts_   = 0;
    std::vector<SCStats> sc_stats_;
    std::vector<int> cluster_sizes_;  // per-output-cluster transcript counts

public:
    MakeClusters(std::vector<ReadList *> &readLists,
                 TranscriptList *tList,
                 std::map<float, std::string> &thresholds,
                 std::vector<int> &groups);

    void set_reading_time(double t) { time_reading_ = t; }
    void print_summary() const;
};