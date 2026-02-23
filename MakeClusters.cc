// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.
//
// Last modified 23 February 2026, Martin Paliocha, martin.paliocha@nmbu.no

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <MakeClusters.h>
#include <Progress.h>

#ifdef HAVE_IGRAPH
#include <LeidenCluster.h>
#endif

using std::cout;
using std::endl;
using std::map;
using std::ofstream;
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

void MakeClusters::makeSuperClusters(vector<ReadList *> &readLists,
                                     TranscriptList *tList) {
    progress::print_banner("Building super-clusters");

    if (!readLists.empty())
        transMap.reserve(1 << 20);  // 1M buckets — will grow if needed

    double t_start = omp_get_wtime();
    int64_t ec_count = 0;
    progress::ProgressLine ec_progress("Super-cluster construction");

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

            ec_count++;
            ec_progress.update(ec_count);
        }
        delete readLists[sample];
    }
    ec_progress.finish();

    readLists.clear();

    // Create singleton clusters for transcripts with only direct counts
    // (no linked reads, so not discovered during union-find above)
    for (auto &[name, transcript] : *tList) {
        if (transcript->has_direct_counts()
                && transMap.find(transcript) == transMap.end()) {
            Cluster *c = new Cluster();
            c->add_tran(transcript);
            c->clusterList_index = static_cast<int>(clusterList.size());
            clusterList.push_back(c);
        }
    }

    transMap.clear();

    time_superclusters_ = omp_get_wtime() - t_start;

    // Sum actual transcript count across all super-clusters
    total_transcripts_ = 0;
    for (Cluster *c : clusterList)
        total_transcripts_ += c->n_trans();

    cout << progress::ansi::green(
            "  " + progress::format_count(static_cast<int64_t>(clusterList.size()))
            + " super-clusters from "
            + progress::format_count(ec_count) + " equivalence classes in "
            + progress::format_duration(time_superclusters_))
         << endl;
}


// ── Stage 2: Parallel clustering ─────────────────────────────────────

void MakeClusters::processSuperClusters(map<float, string> &thresholds,
                                        vector<int> &groups) {
    int n = static_cast<int>(clusterList.size());
    const auto algo = Cluster::algorithm;

    const char *algo_name = "hierarchical";
    if (algo == ClusterAlgorithm::Leiden) algo_name = "Leiden";
    else if (algo == ClusterAlgorithm::Both) algo_name = "hierarchical + Leiden";

    progress::print_banner(string(algo_name) + " clustering ("
                           + progress::format_count(n) + " super-clusters)");

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

    double t_start = omp_get_wtime();
    double last_progress_time = t_start;
    int done = 0;

    // Thread-local SC stats (one vector per thread, merged after loop)
    int max_threads = omp_get_max_threads();
    vector<vector<SCStats>> tl_stats(max_threads);
    vector<vector<int>> tl_csizes(max_threads);

    // Thread-local file buffers: filename → accumulated content.
    // Replaces 750K+ NFS open-append-close operations with 2-4 writes.
    using FileBuffer = map<string, string>;
    vector<FileBuffer> tl_files(max_threads);

    #pragma omp parallel for schedule(dynamic) shared(done, last_progress_time)
    for (int i = 0; i < n; i++) {
        double sc_start = omp_get_wtime();
        Cluster *c = clusterList[i];
        int ntrans = c->n_trans();
        c->set_id(i);
        c->set_sample_groups(groups);

        int nedges = 0;
        vector<int> sizes;
        int tid = omp_get_thread_num();

#ifdef HAVE_IGRAPH
        if (algo == ClusterAlgorithm::Both) {
            // Leiden first (reads cluster data without mutation)
            nedges = cluster_leiden(c, thresholds, "l-", &sizes, &tl_files[tid]);
            // Hierarchical second (mutates via merges — must go last)
            c->cluster(thresholds, "h-", &tl_files[tid]);
        } else if (algo == ClusterAlgorithm::Leiden) {
            nedges = cluster_leiden(c, thresholds, "", &sizes, &tl_files[tid]);
        } else
#endif
        {
            c->cluster(thresholds, "", &tl_files[tid]);
            for (const auto &g : c->get_groups())
                if (!g.empty()) sizes.push_back(static_cast<int>(g.size()));
        }

        double sc_duration = omp_get_wtime() - sc_start;
        tl_stats[tid].push_back({i, ntrans, nedges,
                                 static_cast<int>(sizes.size()), sc_duration});
        auto &cs = tl_csizes[tid];
        cs.insert(cs.end(), sizes.begin(), sizes.end());

        delete c;

        #pragma omp atomic
        done++;

        // Time-based progress: 2s on TTY, 10s on non-TTY
        {
            double now = omp_get_wtime();
            double interval = progress::is_tty() ? 2.0 : 10.0;
            bool should_print = false;
            int done_snapshot = 0;

            #pragma omp critical(progress_time)
            {
                if (now - last_progress_time >= interval) {
                    last_progress_time = now;
                    done_snapshot = done;  // safe read under critical
                    should_print = true;
                }
            }

            if (should_print) {
                double elapsed = now - t_start;
                double rate = (elapsed > 1e-9)
                    ? static_cast<double>(done_snapshot) / elapsed : 0;
                double eta = (rate > 0) ? (n - done_snapshot) / rate : -1;
                double pct = 100.0 * done_snapshot / n;

                char buf[128];
                std::snprintf(buf, sizeof(buf),
                              "[Clustering] %d/%d (%.1f%%) \xe2\x94\x80 %s \xe2\x94\x80 %s remaining",
                              done_snapshot, n, pct,
                              progress::format_rate(rate).c_str(),
                              progress::format_eta(eta).c_str());

                #pragma omp critical(print)
                cout << buf << endl;
            }
        }
    }

    // Merge thread-local stats
    for (auto &tv : tl_stats)
        sc_stats_.insert(sc_stats_.end(), tv.begin(), tv.end());
    for (auto &tv : tl_csizes)
        cluster_sizes_.insert(cluster_sizes_.end(), tv.begin(), tv.end());

    // Merge thread-local file buffers and write to disk in one pass
    {
        FileBuffer merged;
        for (auto &fb : tl_files)
            for (auto &[fn, content] : fb)
                merged[fn] += std::move(content);
        for (auto &[fn, content] : merged) {
            ofstream f(fn);
            f << content;
        }
    }

    time_clustering_ = omp_get_wtime() - t_start;
    clusterList.clear();

    cout << progress::ansi::green(
            "  Clustering complete in " + progress::format_duration(time_clustering_))
         << endl;
}


// ── Summary report ──────────────────────────────────────────────────

// Helper: compute median of a sorted vector
static double median_of(vector<double> &v) {
    if (v.empty()) return 0;
    std::ranges::sort(v);
    size_t m = v.size() / 2;
    return (v.size() % 2 == 0) ? (v[m - 1] + v[m]) / 2.0 : v[m];
}

void MakeClusters::print_summary() const {
    using namespace progress;

    print_banner("Summary");

    // Algorithm description
    const auto algo = Cluster::algorithm;
    string algo_desc;
    if (algo == ClusterAlgorithm::Both)         algo_desc = "hierarchical + Leiden";
    else if (algo == ClusterAlgorithm::Leiden)   algo_desc = "Leiden";
    else                                         algo_desc = "hierarchical";

    // Append Leiden-specific flags
    if (algo != ClusterAlgorithm::Hierarchical) {
        vector<string> flags;
        if (Cluster::knn == 0)       flags.emplace_back("kNN auto");
        else if (Cluster::knn > 0)   flags.push_back("kNN " + std::to_string(Cluster::knn));
        if (Cluster::lrt_softness > 0.0f) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "lrt-softness=%.1f", Cluster::lrt_softness);
            flags.emplace_back(buf);
        }
        if (!flags.empty()) {
            algo_desc += " (";
            for (size_t i = 0; i < flags.size(); i++) {
                if (i > 0) algo_desc += ", ";
                algo_desc += flags[i];
            }
            algo_desc += ")";
        }
    }

    double total = time_reading_ + time_superclusters_ + time_clustering_;

    int64_t total_clusters = static_cast<int64_t>(cluster_sizes_.size());

    std::printf("  Transcripts:     %s\n",
                format_count(total_transcripts_).c_str());
    std::printf("  Super-clusters:  %s\n",
                format_count(static_cast<int64_t>(sc_stats_.size())).c_str());
    std::printf("  Output clusters: %s\n",
                format_count(total_clusters).c_str());
    std::printf("  Algorithm:       %s\n\n", algo_desc.c_str());

    std::printf("  Stage timing:\n");
    std::printf("    Reading input:           %s\n",
                format_duration(time_reading_).c_str());
    std::printf("    Super-cluster formation: %s\n",
                format_duration(time_superclusters_).c_str());
    std::printf("    Clustering:              %s\n",
                format_duration(time_clustering_).c_str());
    std::printf("    Total:                   %s\n", format_duration(total).c_str());

    // ── Output cluster size distribution ──
    if (total_clusters > 0) {
        int64_t n_single = 0, n_small = 0, n_med = 0, n_large = 0, n_xl = 0;
        int max_size = 0;
        for (int sz : cluster_sizes_) {
            if (sz == 1) n_single++;
            else if (sz <= 10) n_small++;
            else if (sz <= 100) n_med++;
            else if (sz <= 1000) n_large++;
            else n_xl++;
            if (sz > max_size) max_size = sz;
        }
        auto sorted_sizes = cluster_sizes_;
        std::ranges::sort(sorted_sizes);
        int median_size = sorted_sizes[sorted_sizes.size() / 2];

        std::printf("\n  Output cluster size distribution:\n");
        std::printf("    %-14s %9s %7s\n", "Size", "Count", "");

        auto print_bucket = [&](const char *label, int64_t count) {
            if (count == 0) return;
            double pct = 100.0 * count / total_clusters;
            std::printf("    %-14s %9s %6.1f%%\n",
                        label, format_count(count).c_str(), pct);
        };
        print_bucket("Singletons", n_single);
        print_bucket("2-10", n_small);
        print_bucket("11-100", n_med);
        print_bucket("101-1K", n_large);
        print_bucket(">1K", n_xl);

        std::printf("    Median: %d   Largest: %s transcripts\n",
                    median_size, format_count(max_size).c_str());
    }

    if (sc_stats_.empty()) {
        std::printf("\n");
        print_separator();
        std::printf("\n");
        std::fflush(stdout);
        return;
    }

    // ── Size-distribution table ──
    struct Bucket {
        const char *label;
        int lo, hi;               // [lo, hi)
        vector<double> durations;
        vector<int64_t> edges;
    };
    Bucket buckets[] = {
        {"1-10",     1,    11,  {}, {}},
        {"11-100",   11,   101, {}, {}},
        {"101-1K",   101,  1001,{}, {}},
        {">1K",      1001, INT_MAX, {}, {}},
    };

    for (const auto &s : sc_stats_) {
        for (auto &b : buckets) {
            if (s.ntrans >= b.lo && s.ntrans < b.hi) {
                b.durations.push_back(s.duration);
                b.edges.push_back(s.nedges);
                break;
            }
        }
    }

    bool has_edges = std::any_of(sc_stats_.begin(), sc_stats_.end(),
                                 [](const SCStats &s) { return s.nedges > 0; });

    std::printf("\n  Super-cluster size distribution:\n");
    if (has_edges) {
        std::printf("    %-10s %7s %10s %10s %12s %12s\n",
                    "Size", "Count", "Med. time", "Max time", "Med. edges", "Max edges");
    } else {
        std::printf("    %-10s %7s %10s %10s\n",
                    "Size", "Count", "Med. time", "Max time");
    }

    for (auto &b : buckets) {
        if (b.durations.empty()) continue;

        double med_dur = median_of(b.durations);
        double max_dur = *std::ranges::max_element(b.durations);
        int64_t count = static_cast<int64_t>(b.durations.size());

        if (has_edges) {
            vector<double> edge_doubles(b.edges.begin(), b.edges.end());
            double med_edges = median_of(edge_doubles);
            int64_t max_edges = *std::ranges::max_element(b.edges);
            std::printf("    %-10s %7s %10s %10s %12s %12s\n",
                        b.label,
                        format_count(count).c_str(),
                        format_duration(med_dur).c_str(),
                        format_duration(max_dur).c_str(),
                        format_count(static_cast<int64_t>(med_edges)).c_str(),
                        format_count(max_edges).c_str());
        } else {
            std::printf("    %-10s %7s %10s %10s\n",
                        b.label,
                        format_count(count).c_str(),
                        format_duration(med_dur).c_str(),
                        format_duration(max_dur).c_str());
        }
    }

    // ── Top 5 slowest super-clusters ──
    auto sorted = sc_stats_;
    std::ranges::sort(sorted, [](const SCStats &a, const SCStats &b) {
        return a.duration > b.duration;
    });

    int top_n = std::min(5, static_cast<int>(sorted.size()));
    std::printf("\n  Top %d slowest super-clusters:\n", top_n);
    for (int i = 0; i < top_n; i++) {
        std::printf("    SC %-6d %7s transcripts   %s",
                    sorted[i].id,
                    format_count(sorted[i].ntrans).c_str(),
                    format_duration(sorted[i].duration).c_str());
        if (sorted[i].nedges > 0)
            std::printf("   %s edges", format_count(sorted[i].nedges).c_str());
        std::printf("\n");
    }

    // ── Bottleneck detection ──
    vector<double> all_durations;
    all_durations.reserve(sc_stats_.size());
    for (const auto &s : sc_stats_) all_durations.push_back(s.duration);
    double med = median_of(all_durations);

    int bottleneck_count = 0;
    if (med > 0.01) {
        for (const auto &s : sc_stats_)
            if (s.duration > 2.0 * med) bottleneck_count++;
    }

    if (bottleneck_count > 0) {
        char buf[128];
        std::snprintf(buf, sizeof(buf),
                      "\n  %s %d super-cluster%s took >2\xc3\x97 median time (%s)",
                      "\xe2\x9a\xa0", bottleneck_count,
                      bottleneck_count > 1 ? "s" : "",
                      format_duration(med).c_str());
        std::printf("%s\n", ansi::yellow(buf).c_str());
    }

    std::printf("\n");
    print_separator();
    std::printf("\n");
    std::fflush(stdout);
}


// ── Constructor: orchestrate both stages ────────────────────────────

MakeClusters::MakeClusters(vector<ReadList *> &readLists,
                           TranscriptList *tList,
                           map<float, string> &thresholds,
                           vector<int> &groups) {
    // Stage 1: group transcripts sharing at least one read
    makeSuperClusters(readLists, tList);
    // Stage 2: cluster within each super-cluster (hierarchical or Leiden)
    processSuperClusters(thresholds, groups);
}
