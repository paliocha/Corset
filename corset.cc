// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** Main program for corset.  Controls I/O and command-line options.
 **
 ** Four main stages:
 ** 1 - Read alignments from BAM / salmon eq_classes / corset files.
 ** 2 - Filter transcripts with fewer than min_counts reads.
 ** 3 - Build super-clusters (transcripts sharing at least one read).
 ** 4 - Clustering within each super-cluster (hierarchical or Leiden, OpenMP parallel).
 **
 ** Original author: Nadia Davidson
 ** OpenMP/htslib port: Martin Paliocha, 2026
 ** Last modified 23 February 2026, martin.paliocha@nmbu.no
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <map>
#include <random>

#include <omp.h>
#include <unistd.h>
#include <htslib/sam.h>

#include <MakeClusters.h>
#include <Progress.h>
#include <Read.h>
#include <Transcript.h>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;

static constexpr std::string_view corset_extension = ".corset-reads";


// ── Read alignments from BAM files (htslib) ─────────────────────────

ReadList *read_bam_file(string all_file_names, TranscriptList *trans, int sample) {
    ReadList *rList = new ReadList(trans);
    string filename;
    stringstream ss(all_file_names);

    while (getline(ss, filename, ',')) {
        cout << "  " << progress::ansi::cyan("[Sample " + std::to_string(sample + 1)
             + "/" + std::to_string(Transcript::samples) + "]")
             << " " << filename << endl;
        if (filename.contains(corset_extension)) {
            cout << "The input files look like corset read files (filenames with "
                 << corset_extension << "). Perhaps you meant to run with the -i corset option?" << endl;
            exit(1);
        }

        samFile *in = sam_open(filename.c_str(), "r");
        if (!in) {
            cerr << "fail to open " << filename << " for reading." << endl;
            exit(1);
        }
        sam_hdr_t *hdr = sam_hdr_read(in);
        if (!hdr) {
            cerr << "fail to read header from " << filename << endl;
            sam_close(in);
            exit(1);
        }
        bam1_t *b = bam_init1();

        // Register all reference sequences from the header
        for (int tid = 0; tid < hdr->n_targets; tid++)
            trans->insert(sam_hdr_tid2name(hdr, tid));

        progress::ProgressLine aln_progress(
            "Sample " + std::to_string(sample + 1) + "/" + std::to_string(Transcript::samples));
        int64_t aln_count = 0;

        while (sam_read1(in, hdr, b) >= 0) {
            string read_name(bam_get_qname(b));
            int tid = b->core.tid;
            if (tid != -1)  // unmapped reads have tid=-1
                rList->add_alignment(read_name, string(sam_hdr_tid2name(hdr, tid)), sample);
            aln_count++;
            aln_progress.update(aln_count);
        }
        aln_progress.finish();

        bam_destroy1(b);
        sam_hdr_destroy(hdr);
        sam_close(in);
    }
    rList->compactify_reads(trans);
    return rList;
}

// ── Read k-mers from FASTA files ────────────────────────────────────

ReadList *read_fasta_file(string all_file_names, TranscriptList *trans, int sample) {
    ReadList *rList = new ReadList(trans);
    string filename;
    stringstream ss(all_file_names);

    while (getline(ss, filename, ',')) {
        cout << "  " << progress::ansi::cyan("[Sample " + std::to_string(sample + 1)
             + "/" + std::to_string(Transcript::samples) + "]")
             << " " << filename << endl;

        ifstream file(filename);
        if (!file.good()) {
            cout << "Unable to open file " << filename << endl;
            exit(1);
        }

        map<string, string> sequences;
        string id, line;
        while (getline(file, line)) {
            int start = static_cast<int>(line.find(">")) + 1;
            if (start == 1) {
                int end = static_cast<int>(line.find_first_of("\t\n ")) - 1;
                id = line.substr(start, end);
            } else {
                sequences[id] += line;
            }
        }

        const int k = 70;
        progress::ProgressLine contig_progress(
            "Sample " + std::to_string(sample + 1) + "/" + std::to_string(Transcript::samples),
            static_cast<int64_t>(sequences.size()));
        int64_t t_count = 0;

        for (auto &kv : sequences) {
            trans->insert(kv.first);
            const string &seq = kv.second;
            for (size_t i = 0; i + k <= seq.length(); i++)
                rList->add_alignment(seq.substr(i, k), kv.first, sample);
            t_count++;
            contig_progress.update(t_count);
        }
        contig_progress.finish();
    }
    rList->compactify_reads(trans);
    return rList;
}


// ── Equivalence class helpers ────────────────────────────────────────

enum class ReadStatus { redistributed, counted, filtered };

// Add an equivalence class into the ReadList.  Common to corset and
// salmon EC file parsing.
ReadStatus add_equivalence_class(ReadList *rList, TranscriptList *tList,
                                 int sample,
                                 vector<string> &transNames, int weight,
                                 std::default_random_engine &rng) {
    if (weight < Transcript::min_reads_for_link) {
        std::shuffle(transNames.begin(), transNames.end(), rng);
        for (size_t rr = 0; rr < transNames.size(); rr++) {
            int this_weight = weight / static_cast<int>(transNames.size());
            if (static_cast<int>(rr) < (weight % static_cast<int>(transNames.size())))
                this_weight++;
            Transcript *t = tList->insert(transNames[rr]);
            t->add_direct_count(sample, this_weight);
        }
        return ReadStatus::redistributed;
    } else if (static_cast<int>(transNames.size()) <= Transcript::max_alignments ||
               Transcript::max_alignments <= 0) {
        rList->add_alignment(transNames, sample, weight);
        return ReadStatus::counted;
    }
    return ReadStatus::filtered;
}


// ── Read corset-format equivalence class files ──────────────────────

ReadList *read_corset_file(string all_file_names, TranscriptList *trans, int sample) {
    ReadList *rList = new ReadList(trans);
    string filename;
    stringstream ss(all_file_names);

    while (getline(ss, filename, ',')) {
        cout << "Reading corset file : " << filename << endl;
        if (!filename.contains(corset_extension)) {
            cout << "The input files don't have the extension, " << corset_extension
                 << ". Please check them." << endl;
            exit(1);
        }

        ifstream file(filename);
        string line;
        int reads_counted = 0, reads_filtered = 0, reads_redistributed = 0;
        std::default_random_engine rng(sample * 2654435761u + 42);

        while (getline(file, line)) {
            istringstream istream(line);
            int weight;
            istream >> weight;
            vector<string> transNames;
            string name;
            while (istream >> name)
                transNames.push_back(name);

            auto ret = add_equivalence_class(rList, trans, sample, transNames, weight, rng);
            switch (ret) {
                case ReadStatus::redistributed: reads_redistributed += weight; break;
                case ReadStatus::counted:       reads_counted += weight;       break;
                case ReadStatus::filtered:      reads_filtered += weight;      break;
            }
        }
        cout << reads_counted << " reads counted, "
             << reads_filtered << " reads filtered, "
             << reads_redistributed << " reads redistributed." << endl;
    }
    return rList;
}

// ── Read salmon equivalence class files ──────────────────────────────
// Pre-resolves transcript name→pointer mapping from the header once,
// then uses O(1) integer indexing instead of string lookups per eq class.

ReadList *read_salmon_eq_classes_file(string all_file_names, TranscriptList *trans, int sample) {
    ReadList *rList = new ReadList(trans);
    string filename;
    stringstream ss(all_file_names);

    while (getline(ss, filename, ',')) {
        cout << "  " << progress::ansi::cyan("[Sample " + std::to_string(sample + 1)
             + "/" + std::to_string(Transcript::samples) + "]")
             << " " << filename << endl;

        ifstream file(filename);
        string line;
        int reads_counted = 0, reads_filtered = 0, reads_redistributed = 0;
        std::default_random_engine rng(sample * 2654435761u + 42);

        // Header: transcript count, then eq class count
        getline(file, line);
        int ntranscripts = std::stoi(line);
        getline(file, line);
        int neqclasses = std::stoi(line);
        cout << "    " << progress::format_count(ntranscripts) << " transcripts, "
             << progress::format_count(neqclasses) << " equivalence classes" << endl;

        // Pre-resolve all transcript names to Transcript* pointers once
        vector<Transcript *> trans_ptrs(ntranscripts);
        for (int nt = 0; nt < ntranscripts; nt++) {
            getline(file, line);
            trans_ptrs[nt] = trans->insert(line);
        }

        // Process equivalence classes
        progress::ProgressLine eq_progress(
            "Sample " + std::to_string(sample + 1) + "/" + std::to_string(Transcript::samples),
            static_cast<int64_t>(neqclasses));

        for (int ne = 0; ne < neqclasses; ne++) {
            getline(file, line);
            istringstream istream(line);
            int eq_size;
            istream >> eq_size;

            vector<Transcript *> eq_trans(eq_size);
            for (int es = 0; es < eq_size; es++) {
                int trans_pos;
                istream >> trans_pos;
                eq_trans[es] = trans_ptrs[trans_pos];
            }
            int weight;
            istream >> weight;

            // Direct counts: redistribute weight without creating Read objects
            if (weight < Transcript::min_reads_for_link) {
                std::shuffle(eq_trans.begin(), eq_trans.end(), rng);
                for (int rr = 0; rr < eq_size; rr++) {
                    int this_weight = weight / eq_size;
                    if (rr < (weight % eq_size))
                        this_weight++;
                    eq_trans[rr]->add_direct_count(sample, this_weight);
                }
                reads_redistributed += weight;
            } else if (eq_size <= Transcript::max_alignments || Transcript::max_alignments <= 0) {
                rList->add_alignment_ptrs(eq_trans, sample, weight);
                reads_counted += weight;
            } else {
                reads_filtered += weight;
            }

            eq_progress.update(ne + 1);
        }
        eq_progress.finish();

        cout << "    " << progress::format_count(reads_counted) << " counted, "
             << progress::format_count(reads_filtered) << " filtered, "
             << progress::format_count(reads_redistributed) << " redistributed" << endl;
    }
    return rList;
}


// ── Usage text ──────────────────────────────────────────────────────

void print_usage() {
    cout << "\n"
         << "Corset clusters contigs and counts reads from de novo assembled transcriptomes.\n"
         << "\n"
         << "Usage: corset [options] <input bam files>\n"
         << "\n"
         << "Input bam files:\n"
         << "  The input files should be multi-mapped bam files. They can be single, paired-end\n"
         << "  or mixed and do not need to be indexed. A space separated list should be given.\n"
         << "  e.g. corset sample1.bam sample2.bam sample3.bam\n"
         << "  or just: corset sample*.bam\n"
         << "\n"
         << "  If you want to combine the results from different transcriptomes (same reads\n"
         << "  mapped twice or more), use a comma separated list:\n"
         << "  corset sample1_Trinity.bam,sample1_Oases.bam sample2_Trinity.bam,sample2_Oases.bam ...\n"
         << "\n"
         << "Options:\n"
         << "\n"
         << "  -d <double list>  Comma-separated distance thresholds (0-1). Default: 0.3\n"
         << "  -D <double>       Log likelihood ratio threshold. Default: 17.5 + 2.5*ndf\n"
         << "  -I                Switch off the log likelihood ratio test.\n"
         << "  -m <int>          Filter transcripts with fewer than this many reads. Default: 10\n"
         << "  -g <list>         Comma-separated sample groupings.\n"
         << "  -p <string>       Output filename prefix.\n"
         << "  -f <true/false>   Overwrite existing output files. Default: false\n"
         << "  -n <string list>  Comma-separated sample names for count file header.\n"
         << "  -r <true/true-stop/false>  Output read alignment summary. Default: false\n"
         << "  -i <bam/corset/salmon_eq_classes>  Input file type. Default: bam\n"
         << "  -l <int>          Min reads for a link (corset/salmon_eq_classes mode). Default: 1\n"
         << "  -x <int>          Max alignments per read (corset/salmon_eq_classes mode).\n"
         << "  -t <int>          Threads for parallel clustering. Default: auto\n"
         << "  --algorithm <str> Clustering algorithm: hierarchical (default), leiden, or both.\n"
         << "                    Leiden uses -d thresholds as CPM resolution values.\n"
         << "                    'both' outputs h- and l- prefixed files for comparison.\n"
         << "  --lrt-softness <float>  Sigmoid steepness for soft LRT (Leiden only). Default: 0 (hard cutoff).\n"
         << "                    Recommended starting point: 2.0 for gentle sigmoid decay.\n"
         << "  --knn <int|auto>  kNN graph sparsification (Leiden only). Default: auto.\n"
         << "  --no-knn          Disable kNN sparsification (use full edge graph).\n"
         << "                    'auto' = quality-preserving adaptive k per super-cluster.\n"
         << "                    <int> = fixed global k for all super-clusters.\n"
         << "  -v, --version     Print version and exit.\n"
         << "  -h, --help        Print this help message and exit.\n"
         << "\n"
         << "Citation: Nadia M. Davidson and Alicia Oshlack, Corset: enabling differential gene\n"
         << "          expression analysis for de novo assembled transcriptomes, Genome Biology 2014, 15:410\n"
         << "\n"
         << "Please see https://github.com/Oshlack/Corset/wiki for more information\n"
         << endl;
}

// ── Main ────────────────────────────────────────────────────────────

int main(int argc, char **argv) {
    // Default parameters
    string distance_string("0.3");
    bool force = false;
    string sample_names;
    int c;
    vector<int> groups;
    bool output_reads = false;
    bool stop_after_read = false;

    // Function pointer to the input reader (BAM by default)
    ReadList *(*read_input)(string, TranscriptList *, int) = read_bam_file;

    // Handle long options before getopt (which only does short opts).
    // Mark consumed args with nullptr, then compact argv.
    for (int i = 1; i < argc; ++i) {
        string arg(argv[i]);
        if (arg == "--version" || arg == "-v") {
            cout << "Corset version " << VERSION << endl;
            return 0;
        }
        if (arg == "--help" || arg == "-h") {
            print_usage();
            return 0;
        }
        if (arg == "--algorithm" && i + 1 < argc) {
            string val(argv[i + 1]);
            if (val == "hierarchical") {
                Cluster::algorithm = ClusterAlgorithm::Hierarchical;
            } else if (val == "leiden" || val == "both") {
#ifdef HAVE_IGRAPH
                Cluster::algorithm = (val == "leiden")
                    ? ClusterAlgorithm::Leiden
                    : ClusterAlgorithm::Both;
#else
                cerr << "ERROR: " << val << " requires igraph. Rebuild with --with-igraph=DIR." << endl;
                exit(1);
#endif
            } else {
                cerr << "ERROR: Unknown algorithm '" << val
                     << "'. Choose 'hierarchical', 'leiden', or 'both'." << endl;
                exit(1);
            }
            argv[i] = nullptr; argv[++i] = nullptr;
            continue;
        }
        if (arg == "--lrt-softness" && i + 1 < argc) {
            Cluster::lrt_softness = std::stof(string(argv[i + 1]));
            argv[i] = nullptr; argv[++i] = nullptr;
            continue;
        }
        if (arg == "--knn" && i + 1 < argc) {
            string val(argv[i + 1]);
            if (val == "auto")
                Cluster::knn = 0;
            else
                Cluster::knn = std::stoi(val);
            argv[i] = nullptr; argv[++i] = nullptr;
            continue;
        }
        if (arg == "--no-knn") {
            Cluster::knn = -1;
            argv[i] = nullptr;
            continue;
        }
    }
    // Compact argv: remove consumed long options
    {
        int dst = 1;
        for (int src = 1; src < argc; src++)
            if (argv[src]) argv[dst++] = argv[src];
        argc = dst;
    }

    cout << "\n" << progress::ansi::bold("Corset " + string(VERSION)) << endl;
    cout << "  Threads:    " << omp_get_max_threads() << endl;
    if (Cluster::algorithm == ClusterAlgorithm::Both)
        cout << "  Algorithm:  both (hierarchical + Leiden)" << endl;
    else if (Cluster::algorithm == ClusterAlgorithm::Leiden)
        cout << "  Algorithm:  Leiden (CPM)" << endl;
    else
        cout << "  Algorithm:  hierarchical" << endl;

    // Warn if Leiden-only flags used with hierarchical-only mode
    if (Cluster::algorithm == ClusterAlgorithm::Hierarchical) {
        if (Cluster::lrt_softness > 0.0f) {
            cerr << "WARNING: --lrt-softness ignored (only applies to Leiden)" << endl;
            Cluster::lrt_softness = 0.0f;
        }
        Cluster::knn = -1;  // kNN is a no-op for hierarchical
    }

    if (Cluster::lrt_softness > 0.0f)
        cout << "  LRT softness: " << Cluster::lrt_softness << " (sigmoid decay)" << endl;
    if (Cluster::knn == 0)
        cout << "  kNN: auto (quality-preserving per super-cluster)" << endl;
    else if (Cluster::knn > 0)
        cout << "  kNN: " << Cluster::knn << " (fixed)" << endl;

    // Allow two levels of parallelism: outer (super-clusters) + inner (merge)
    omp_set_max_active_levels(2);

    // Parse command-line options
    while ((c = getopt(argc, argv, "f:p:d:n:g:D:Im:r:i:l:x:t:")) != EOF) {
        switch (c) {
        case 'f': {
            string value(optarg);
            std::ranges::transform(value, value.begin(), ::tolower);
            if (value == "true" || value == "t" || value == "1") {
                force = true;
                cout << "Setting output files to be overridden" << endl;
            } else if (value == "false" || value == "f" || value == "0") {
                force = false;
            } else {
                cerr << "Unknown argument passed with -f. Please specify true or false." << endl;
                print_usage();
                exit(1);
            }
            break;
        }
        case 'p':
            cout << "Setting output filename prefix to " << optarg << endl;
            Cluster::file_prefix = string(optarg) + "-";
            break;
        case 'd':
            cout << "Setting distance threshold to: " << optarg << endl;
            distance_string = optarg;
            break;
        case 'n': {
            cout << "Setting sample names to:" << optarg << endl;
            sample_names = string(optarg);
            std::ranges::replace(sample_names, ',', '\t');
            break;
        }
        case 'g': {
            cout << "Setting sample groups:" << optarg;
            stringstream ss(optarg);
            string s;
            vector<string> names;
            while (getline(ss, s, ','))
                names.push_back(s);
            int ngroups = 0;
            for (size_t s1 = 0; s1 < names.size(); s1++) {
                size_t s2 = 0;
                for (; s2 < names.size(); s2++) {
                    if (names[s1] == names[s2])
                        break;
                }
                if (groups.size() <= s2) {
                    groups.push_back(ngroups);
                    ngroups++;
                } else {
                    groups.push_back(groups[s2]);
                }
            }
            Transcript::groups = ngroups;
            cout << ", " << Transcript::groups << " groups in total" << endl;
            break;
        }
        case 'D':
            cout << "Setting likelihood threshold at " << optarg << endl;
            Cluster::D_cut = std::stof(string(optarg));
            break;
        case 'I':
            cout << "Switching likelihood test off" << endl;
            Cluster::D_cut = INT_MAX;
            break;
        case 'm':
            cout << "Setting minimum counts to " << optarg << endl;
            Transcript::min_counts = std::stoi(string(optarg));
            break;
        case 'r': {
            string value(optarg);
            std::ranges::transform(value, value.begin(), ::tolower);
            if (value == "true-stop") {
                output_reads = true;
                stop_after_read = true;
                cout << "Setting read alignments to be output to file.\n"
                     << "Corset will exit after reading bam files." << endl;
            } else if (value == "true" || value == "t") {
                output_reads = true;
                cout << "Setting read alignments to be output to file." << endl;
            } else if (value != "false" && value != "f") {
                cerr << "Unknown argument passed with -r. Please specify true or false." << endl;
                print_usage();
                exit(1);
            }
            break;
        }
        case 'i': {
            string value(optarg);
            std::ranges::transform(value, value.begin(), ::tolower);
            if (value == "corset") {
                read_input = read_corset_file;
                output_reads = false;
                stop_after_read = false;
            } else if (value == "salmon_eq_classes") {
                read_input = read_salmon_eq_classes_file;
                output_reads = false;
                stop_after_read = false;
            } else if (value == "fasta") {
                read_input = read_fasta_file;
            } else if (value != "bam") {
                cerr << "Unknown input type, " << value
                     << ", passed with -i. Please check options." << endl;
                print_usage();
                exit(1);
            }
            break;
        }
        case 'l':
            cout << "Setting minimum reads for a link to " << optarg << endl;
            Transcript::min_reads_for_link = std::stoi(string(optarg));
            break;
        case 'x':
            cout << "Setting maximum alignments for a read to " << optarg << endl;
            Transcript::max_alignments = std::stoi(string(optarg));
            break;
        case 't': {
            int nthreads = std::stoi(string(optarg));
            if (nthreads < 1) nthreads = 1;
            omp_set_num_threads(nthreads);
            cout << "Overriding thread count to " << nthreads << endl;
            break;
        }
        case '?':
            cerr << "Unknown option.. stopping" << endl;
            print_usage();
            exit(1);
            break;
        }
    }

    // Parse distance thresholds
    map<float, string> distance_thresholds;
    {
        istringstream dss(distance_string);
        string token;
        while (getline(dss, token, ',')) {
            float value = std::stof(token);
            if (value < 0 || value > 1) {
                cerr << "The distance provided is invalid. Must be between 0 and 1." << endl;
                print_usage();
                exit(1);
            }
            distance_thresholds[value] = string("-") + token;
        }
    }
    if (distance_thresholds.size() == 1)
        distance_thresholds.begin()->second = "";

    // Remaining arguments are input files
    int smpls = argc - optind;
    if (smpls == 0) {
        cerr << "No input files specified" << endl;
        print_usage();
        exit(1);
    }
    if (smpls > 255) {
        cerr << "Error: maximum of 255 samples supported (got " << smpls << ")." << endl;
        exit(1);
    }

    if (sample_names.empty()) {
        for (int s = optind; s < argc; s++) {
            sample_names += string(argv[s]);
            if (s < argc - 1) sample_names += '\t';
        }
    }

    int n_sample_names = static_cast<int>(std::ranges::count(sample_names, '\t')) + 1;
    if (n_sample_names != smpls) {
        cerr << "The number of sample names passed (via -n), " << n_sample_names
             << ", does not match the number of samples, " << smpls << "." << endl;
        exit(1);
    }

    if (groups.empty()) {
        for (int s = 0; s < smpls; s++)
            groups.push_back(s);
        Transcript::groups = smpls;
    }

    if (static_cast<int>(groups.size()) != smpls) {
        cerr << "The number of experimental groups (via -g) does not match "
             << "the number of input files." << endl;
        exit(1);
    }

    // Default D_cut based on degrees of freedom
    if (Cluster::D_cut == 0) {
        Cluster::D_cut = 17.5f + 2.5f * (Transcript::groups - 1);
        cout << "  D threshold:  " << Cluster::D_cut
             << " (default for " << Transcript::groups << " groups)" << endl;
    }

    // Determine method tags for output files
    vector<string> method_tags;
    if (Cluster::algorithm == ClusterAlgorithm::Both)
        method_tags = {"h-", "l-"};
    else
        method_tags = {""};

    // Check / create output files
    for (auto it = distance_thresholds.begin();
         !stop_after_read && it != distance_thresholds.end(); ++it) {
        std::string_view types[] = {Cluster::file_counts, Cluster::file_clusters};
        for (const auto &tag : method_tags) {
            for (int t = 0; t < 2; t++) {
                string fn = Cluster::file_prefix + tag + string(types[t])
                          + it->second + string(Cluster::file_ext);
                ifstream probe(fn);
                if (probe.good()) {
                    if (force && remove(fn.c_str()) != 0) {
                        cerr << "Could not replace the file, " << fn << endl;
                        exit(1);
                    } else if (!force) {
                        cerr << "File already exists, " << fn
                             << ". Use \"-f true\" to overwrite." << endl;
                        exit(1);
                    }
                }
                probe.close();
                ofstream ofile(fn);
                if (types[t] == Cluster::file_counts)
                    ofile << '\t' + sample_names << endl;
                ofile.close();
            }
        }
    }

    // Read input files
    Transcript::samples = smpls;
    TranscriptList *tList = new TranscriptList;
    vector<ReadList *> rList;

    progress::print_banner("Reading input files (" + std::to_string(smpls) + " samples)");
    double t_read_start = omp_get_wtime();

    for (int f = 0; f < smpls; f++) {
        rList.push_back(read_input(string(argv[optind + f]), tList, f));
        if (output_reads) {
            string bam_fn = string(argv[optind + f]);
            string base_fn = bam_fn.substr(bam_fn.find_last_of("/\\,") + 1);
            rList[f]->write(base_fn + string(corset_extension));
        }
    }

    double reading_time = omp_get_wtime() - t_read_start;
    cout << progress::ansi::green(
            "  All " + std::to_string(smpls) + " samples read in "
            + progress::format_duration(reading_time))
         << endl;
    if (stop_after_read) exit(0);

    // Handle transcripts below minimum counts
    int n = 0;
    for (auto &[name, transcript] : *tList) {
        if (!transcript->reached_min_counts()) {
            if (Transcript::min_counts == 0) {
                n++;
                for (auto &dis : distance_thresholds) {
                    for (const auto &tag : method_tags) {
                        string fn = Cluster::file_prefix + tag
                                    + string(Cluster::file_clusters)
                                    + dis.second + string(Cluster::file_ext);
                        ofstream clusterFile(fn, std::ios_base::app);
                        clusterFile << transcript->get_name() << "\t"
                                    << Cluster::cluster_id_prefix_no_reads << n << "\n";
                    }
                }
            }
            transcript->remove();
        }
    }

    MakeClusters cList(rList, tList, distance_thresholds, groups);
    cList.set_reading_time(reading_time);
    cList.print_summary();

    delete tList;
    return 0;
}