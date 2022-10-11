#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

#define LENGTH 1000000000
#define KMER_LENGTH 31

using namespace std;

void print_usage();

void create_kmer(char kmer[], int length, const string& nuc_set);

int main(int argc, char** argv ) {
    cout << "### Random sequence generator ###" << endl;
    // parameters
    int seed = time(nullptr);
    int length = LENGTH;
    int kmer_length = KMER_LENGTH;

    // get parameters if any
    for(int i = 0; i < argc; i++){
        if((strcmp(argv[i], "-k") == 0) && (i < argc - 1))
            kmer_length = stoi(argv[++i]);
        if((strcmp(argv[i], "-s") == 0) && (i < argc - 1))
            seed = stoi(argv[++i]);
        if((strcmp(argv[i], "-l") == 0) && (i < argc - 1))
            length = stoi(argv[++i]);
        if((strcmp(argv[i], "-h") == 0)) {
            print_usage();
            exit(0);
        }
    }

    // print acquired parameters
    cout << "Parameters: seed=" << seed << "; length=" << length << "; kmer_length=" << kmer_length << endl;

    srand(seed);
    int num_kmers = length / kmer_length;
    string filename = "RND"+to_string(seed)+".fasta";
    fstream fout;
    fout.open(filename, ios::out);

    cout << setw(5) << setprecision(3);
    char buf[kmer_length + 3];
    buf[0] = '>';
    buf[1] = '\n';
    buf[kmer_length + 2] = '\n';
    for(int line = 0; line < num_kmers; line++){
        create_kmer(buf + 2, kmer_length, "ACTG");
        fout.write(buf, kmer_length + 3);

        // progress bar
        if(line % 1000000 == 0) {
            cout << "writing " << 1.0 * line / num_kmers * 100<< " %...\r";
            cout.flush();
        }
    }
    cout << "writing   100%" << endl;
    fout.close();
    cout << "Wrote " << num_kmers << " random kmers to " << filename << endl;
    return 0;
}

void print_usage() {
    cout    << "Create a FASTA file with random reads.\nParameters:\n";
    cout    << "\t-k\t\tkmer-length [" << KMER_LENGTH <<"]\n"
            << "\t-l\t\tfile size [" << LENGTH << "]\n"
            << "\t-s\t\trandom seed [seconds since Epoch]" << endl;
}

void create_kmer(char kmer[], int length, const string& nuc_set) {
    for(int i = 0; i < length; i++){
        int n = rand() % nuc_set.length();
        // equiprobable nucleotides with nuc_set="ACTG"
        kmer[i] = nuc_set[n];
    }
}
