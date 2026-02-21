// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// Builds super-clusters (stage 1) via weighted union-find, then
// dispatches parallel hierarchical clustering (stage 2).
//
// Original author: Nadia Davidson
// Last modified 21 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <Cluster.h>
#include <Read.h>
#include <Transcript.h>

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

    void makeSuperClusters(std::vector<ReadList *> &readLists);
    void processSuperClusters(std::map<float, std::string> &thresholds,
                              std::vector<int> &groups);

public:
    MakeClusters(std::vector<ReadList *> &readLists,
                 std::map<float, std::string> &thresholds,
                 std::vector<int> &groups);
};