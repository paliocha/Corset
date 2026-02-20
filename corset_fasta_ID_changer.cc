// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

// Prefix FASTA contig IDs with their cluster ID from a corset
// clusters file.  Thanks to Marco Salvemini for suggesting this.
//
// Usage: corset_fasta_ID_changer <cluster file> <fasta file> > <out>
//
// Author: Nadia Davidson

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cstdlib>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::map;
using std::string;

static void print_usage() {
    cout << "\ncorset_fasta_ID_changer modifies a fasta file by prefixing "
         << "each contig ID with its associated cluster ID\n\n"
         << "Usage: corset_fasta_ID_changer <cluster file> <fasta file>  >  <out file>\n"
         << "  <cluster file> - Table with two columns: contig ID and cluster ID.\n"
         << "  <fasta file>   - FASTA file whose contig IDs should be altered.\n"
         << "  <out file>     - Output FASTA to stdout.\n\n"
         << "See https://github.com/Oshlack/Corset/wiki for more information\n"
         << endl;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        print_usage();
        exit(1);
    }

    // Read cluster→contig mapping
    map<string, string> cluster_contig_map;
    ifstream file(argv[1]);
    if (!file.good()) {
        cerr << "Unable to open file " << argv[1] << endl;
        exit(1);
    }

    string line, contig, cluster;
    while (getline(file, line)) {
        istringstream ls(line);
        ls >> contig >> cluster;
        cluster_contig_map[contig] = cluster;
    }
    file.close();

    // Process FASTA, prefixing IDs
    ifstream fasta(argv[2]);
    if (!fasta.good()) {
        cerr << "Unable to open file " << argv[2] << endl;
        exit(1);
    }

    while (getline(fasta, line)) {
        if (line.starts_with('>')) {
            int end = static_cast<int>(line.find_first_of("\t\n ")) - 1;
            string id = line.substr(1, end);
            string cluster_id = cluster_contig_map[id];
            if (cluster_id.empty())
                cluster_id = "UnknownCluster";
            line.replace(1, end, cluster_id + "--" + id);
        }
        cout << line << "\n";
    }
    fasta.close();

    return 0;
}
