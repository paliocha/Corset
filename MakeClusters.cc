// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.
//
// Last modified 22 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <MakeClusters.h>

#ifdef HAVE_IGRAPH
#include <LeidenCluster.h>
#endif

using std::cout;
using std::endl;
using std::map;
using std::string;
using std::vector;


// ── Find-or-create the cluster for a transcript ─────────────────────

std::pair<Transcript *const, Cluster *> *
MakeClusters::getMapElement(Transcript *trans) {
    auto it = transMap.find(trans);
    if (it != transMap.end())
        return &(*it);

    Cluster *clust = new Cluster();
    clust->add_tran(trans);
    clust->clusterList_index = static_cast<int>(clusterList.size());
    clusterList.push_back(clust);

    auto result = transMap.emplace(trans, clust);
    return &(*result.first);
}


// ── Weighted union-find merge ───────────────────────────────────────
// Always absorb smaller cluster into larger to guarantee each
// transcript changes cluster at most O(log n) times total.

void MakeClusters::checkAgainstCurrentCluster(Transcript *trans) {
    Cluster *this_cluster = getMapElement(trans)->second;
    if (this_cluster == current_cluster) return;

    Cluster *big   = current_cluster;
    Cluster *small = this_cluster;
    if (small->n_trans() > big->n_trans()) {
        std::swap(big, small);
        current_cluster = big;
    }

    // Move all transcripts from small → big
    for (int i = 0; i < small->n_trans(); i++) {
        transMap[small->get_tran(i)] = big;
        big->add_tran(small->get_tran(i));
    }
    // Copy all reads across
    for (int i = 0; i < small->n_reads(); i++)
        big->add_read(small->get_read(i));

    // O(1) removal: swap with back and pop
    int idx = small->clusterList_index;
    Cluster *back = clusterList.back();
    clusterList[idx] = back;
    back->clusterList_index = idx;
    clusterList.pop_back();
    delete small;
}


// ── Stage 1: Build super-clusters from shared reads ─────────────────

void MakeClusters::makeSuperClusters(vector<ReadList *> &readLists) {
    if (!readLists.empty())
        transMap.reserve(1 << 20);  // 1M buckets — will grow if needed

    int ec_count = 0;
    for (size_t sample = 0; sample < readLists.size(); sample++) {
        ReadList reads = *(readLists[sample]);

        for (Read *r : reads) {
            auto tIt     = r->align_begin();
            auto tIt_end = r->align_end();

            for (auto t = tIt; t != tIt_end; ++t) {
                if (t == tIt) {
                    setCurrentCluster(*t);
                    current_cluster->add_read(r);
                } else {
                    checkAgainstCurrentCluster(*t);
                }
            }

            if (ec_count % 100000 == 0)
                cout << static_cast<float>(ec_count) / 1e6f
                     << " million equivalence classes read" << endl;
            ec_count++;
        }
        delete readLists[sample];
    }

    readLists.clear();
    transMap.clear();
}


// ── Stage 2: Parallel clustering ─────────────────────────────────────

void MakeClusters::processSuperClusters(map<float, string> &thresholds,
                                        vector<int> &groups) {
    int n = static_cast<int>(clusterList.size());
    const auto algo = Cluster::algorithm;

    const char *algo_name = "hierarchical";
    if (algo == ClusterAlgorithm::Leiden) algo_name = "Leiden";
    else if (algo == ClusterAlgorithm::Both) algo_name = "hierarchical + Leiden (comparison)";
    cout << "Starting " << algo_name << " clustering of "
         << n << " super-clusters..." << endl;

    // Sort largest-first so OMP dynamic schedule starts monster clusters
    // early, preventing a long tail where one thread processes the last
    // giant super-cluster while all others are idle.
    // Only for Leiden/Both: hierarchical allocates O(n^2) per-cluster
    // memory, so launching the 32 largest simultaneously causes swapping.
    if (algo != ClusterAlgorithm::Hierarchical) {
        std::ranges::sort(clusterList, [](Cluster *a, Cluster *b) {
            return a->n_trans() > b->n_trans();
        });
    }

    int done = 0;
    #pragma omp parallel for schedule(dynamic) shared(done)
    for (int i = 0; i < n; i++) {
        Cluster *c = clusterList[i];
        c->set_id(i);
        c->set_sample_groups(groups);

#ifdef HAVE_IGRAPH
        if (algo == ClusterAlgorithm::Both) {
            // Leiden first (reads cluster data without mutation)
            cluster_leiden(c, thresholds, "l-");
            // Hierarchical second (mutates via merges — must go last)
            c->cluster(thresholds, "h-");
        } else if (algo == ClusterAlgorithm::Leiden) {
            cluster_leiden(c, thresholds);
        } else
#endif
        {
            c->cluster(thresholds);
        }

        delete c;

        #pragma omp atomic
        done++;
        if (done % 1000 == 0) {
            #pragma omp critical(print)
            cout << done << " of " << n << " super-clusters done" << endl;
        }
    }
    clusterList.clear();
}


// ── Constructor: orchestrate both stages ────────────────────────────

MakeClusters::MakeClusters(vector<ReadList *> &readLists,
                           map<float, string> &thresholds,
                           vector<int> &groups) {
    // Stage 1: group transcripts sharing at least one read
    makeSuperClusters(readLists);
    // Stage 2: cluster within each super-cluster (hierarchical or Leiden)
    processSuperClusters(thresholds, groups);
}
