#!/bin/env python3
import sys

import pandas as pd
import matplotlib.pyplot as plt

# DEFAULT_RESULTS = "results/results_2022-10-04_14.41.16.csv"
DEFAULT_RESULTS = "results.csv"
SEQUENCES_FILE = "sequences-test.txt"
# methods: output_list
METHODS = {"none": ["fasta"], "bcalm": ["fasta"], "prophasm": ["fasta"], "ust": ["fasta", "counts"],
           "metagraph": ["dbg", "annodbg", "counts"], "contigs": ["fasta", "kmer_counts"]}
NO_COUNTS_METHOD = ["prophasm"]
FIGURES_PATH = "figures"


def read_from_file(myfile):
    seq = []
    with open(myfile, "r") as file:
        for line in file:
            sequence = line.strip()
            # skip white lines or comments
            if sequence == "" or sequence[0] == "#":
                continue
            seq.append(sequence)
    return seq


# function to add value labels
def addlabels2(x, y1, y2):
    for i in range(len(x)):
        plt.text(i, y1[i], f"{y1[i] * 1.0:.2}", ha='center', va='bottom', color='green')
        plt.text(i, y2[i], f"{y2[i] * 1.0:.2}", ha='center', va='top', color='green')


def addlabels1(x, y):
    for i in range(len(x)):
        plt.text(i, y[i], f"{y[i]:.2}", ha='center', va='bottom', color='green')


def analyze(sequences: list, kmer_sizes: list, data: pd.DataFrame, counts=False):
    if counts:
        counts_name = "counts"
        inv_counts_name = "no-counts"
        methods = [e for e in list(METHODS) if e not in NO_COUNTS_METHOD]
    else:
        counts_name = "no-counts"
        inv_counts_name = "counts"
        methods = list(METHODS)
    seq_sizes = \
        [data.query(f'sequence == "{sequence}" and method == "none" and compression == "none"')['size'].iloc[0]
         for sequence in sequences]
    seq_sizes = dict(zip(sequences, seq_sizes))

    print(counts_name)
    all_sizes = {}
    compression_stat = {}
    for kmer_size in kmer_sizes:
        all_sequences = {}
        print(f"\tkmer-size: {kmer_size}")
        avg_ratios = [0 for i in range(len(methods))]
        avg_ratios_ = [0 for i in range(len(methods))]
        avg_size = [0 for i in range(len(methods))]
        avg_size_ = [0 for i in range(len(methods))]
        # average compressed size
        avg_comp = [0 for i in range(len(methods))]
        for sequence in sequences:
            print(f"\t\tsequence: {sequence}")
            # select sequence and kmer-size
            seq_df = data.query(
                f'sequence == "{sequence}" and counts != "{inv_counts_name}" and kmer_size == {kmer_size}'
            )

            # for each method find
            # - the size
            # - best compressed size and tool
            no_compression_sizes = []
            best_compression_sizes = []
            # best_compression_tool = []
            for i, method in enumerate(methods):
                # sum all sizes of uncompressed files for this method
                size = seq_df.query(f'method == "{method}" and compression == "none"')['size'].sum()
                no_compression_sizes.append(size)

                # find the best compression size
                comp_size = 0
                # consider all file types
                for file_type in METHODS.get(method):
                    q = seq_df.query(f'method == "{method}" and compression != "none" and file_type == "{file_type}"')
                    comp_size += q['size'].min()

                    # find the best compression tool
                    idx = q['size'].idxmin()
                    compression = data['compression'].iloc[idx]
                    compression_stat.update({compression: 1 + compression_stat.get(compression, 0)})

                    # only one file if no counts
                    if not counts:
                        break

                best_compression_sizes.append(comp_size)

                avg_ratios[i] += (size / comp_size) / len(sequences)
                avg_ratios_[i] += (seq_sizes[sequence] / comp_size) / len(sequences)
                avg_size[i] += size / len(sequences)
                avg_comp[i] += comp_size / len(sequences)

            # methods results
            print(f"\t\t\tstarting size: {seq_sizes[sequence]}")
            print(f"\t\t\tuncompressed: {no_compression_sizes}")
            print(f"\t\t\tcompressed: {best_compression_sizes}")
            print(
                f"\t\t\tratio1: {[seq_sizes[sequence] / no_compression_size for no_compression_size in no_compression_sizes]}")
            print(
                f"\t\t\tratio2: {[seq_sizes[sequence] / best_compression_size for best_compression_size in best_compression_sizes]}")

            # plot methods vs sizes
            plt.title(f"{sequence} - {counts_name} - k={kmer_size}")
            plt.xlabel("methods")
            plt.ylabel("sizes [bytes]")
            plt.bar(methods, no_compression_sizes)
            plt.bar(methods, best_compression_sizes)
            addlabels2(methods, no_compression_sizes, best_compression_sizes)
            plt.savefig(f"{FIGURES_PATH}/{sequence}.{counts_name}.k{kmer_size}.png")
            plt.clf()

            # save sizes
            all_sequences.update({sequence: [no_compression_sizes, best_compression_sizes]})

        # plot average ratios
        plt.title(f"average - {counts_name} - k={kmer_size}")
        plt.xlabel("methods")
        plt.ylabel("size ratios bis")
        plt.bar(methods, avg_ratios)
        addlabels1(methods, avg_ratios)
        plt.savefig(f"{FIGURES_PATH}/average-ratios-bis.{counts_name}.k{kmer_size}.png")
        plt.clf()

        # plot average ratios
        plt.title(f"average - {counts_name} - k={kmer_size}")
        plt.xlabel("methods")
        plt.ylabel("compression ratios")
        plt.bar(methods, avg_ratios_)
        addlabels1(methods, avg_ratios_)
        plt.savefig(f"{FIGURES_PATH}/average-ratios.{counts_name}.k{kmer_size}.png")
        plt.clf()

        # plot average (compressed) sizes
        plt.title(f"average - {counts_name} - k={kmer_size}")
        plt.xlabel("methods")
        plt.ylabel("sizes [bytes]")
        plt.bar(methods, avg_size)
        plt.bar(methods, avg_comp)
        addlabels2(methods, avg_size, avg_comp)
        plt.savefig(f"{FIGURES_PATH}/average.{counts_name}.k{kmer_size}.png")
        plt.clf()

        all_sizes.update({kmer_size: all_sequences})
    print(compression_stat)
    return all_sizes


def main():
    print("*** Executing plot utility...")
    # parse parameter or set default
    if len(sys.argv) == 2:
        infile = sys.argv[1]
    else:
        infile = DEFAULT_RESULTS

    # read csv file
    # headers: sequence,method,counts,kmer_size,file_type,compression,size
    data = pd.read_csv(infile)

    # read sequence file
    sequences = read_from_file(SEQUENCES_FILE)

    # read kmer-sizes file
    kmer_sizes_file = "kmer-sizes.txt"
    kmer_sizes = read_from_file(kmer_sizes_file)

    # no-counts
    analyze(sequences, kmer_sizes, data, counts=False)
    # counts
    analyze(sequences, kmer_sizes, data, counts=True)


if __name__ == "__main__":
    main()
