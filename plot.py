#!/bin/env python3
import sys

import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_RESULTS = "results/results_2022-09-08_11.42.35.csv"
SEQUENCES_FILE = "sequences-test.txt"
# methods: output_list
METHODS = {"bcalm": ["fasta"], "prophasm": ["fasta"], "ust": ["fasta", "counts"],
           "metagraph": ["dbg", "annodbg", "coords"], "assembly": ["fasta"]}
NO_COUNTS_METHOD = ["prophasm", "assembly"]
FIGURES_PATH = "figures"


def read_from_file(myfile):
    seq = []
    with open(myfile, "r") as file:
        for line in file:
            sequence = line.strip()
            if sequence == "" or sequence[0] == "#":
                continue
            seq.append(sequence)
    return seq


# function to add value labels
def addlabels(x, y1, y2):
    for i in range(len(x)):
        # plt.text(i, y1[i], f"{y1[i]} -> {y2[i]}", ha='center', va='top') # why top does not work?
        plt.text(i, y1[i], f"{y1[i]}", ha='center', va='bottom', color='green')
        plt.text(i, y2[i], f"{y2[i]}", ha='center', va='top', color='green')


def plot_sequences_vs_methods(sequences: list, kmer_sizes: list, data: pd.DataFrame, counts: bool):
    if counts:
        counts_name = "counts"
        methods = [e for e in list(METHODS) if e not in NO_COUNTS_METHOD]
    else:
        counts_name = "no-counts"
        methods = list(METHODS)
    for kmer_size in kmer_sizes:
        print(f"kmer-size: {kmer_size}")
        for sequence in sequences:
            print(f"\tsequence: {sequence}")
            # select sequence and kmer-size
            seq_df = data.query(
                f'sequence == "{sequence}" and counts == "{counts_name}" and kmer_size == {kmer_size}'
            )

            # for each method find
            # - the size
            # - best compressed size and tool
            no_compression_size = []
            best_compression_size = []
            # best_compression_tool = []
            for method in methods:
                # sum all sizes of uncompressed files for this method
                size = seq_df.query(f'method == "{method}" and compression == "none"')['size'].sum()
                no_compression_size.append(size)

                # find the best compression size
                size = 0
                for file_type in METHODS.get(method):
                    q = seq_df.query(f'method == "{method}" and compression != "none" and file_type == "{file_type}"')
                    size += q['size'].min()
                    if not counts:
                        break
                best_compression_size.append(size)

                # find the best compression tool
                # idx = q['size'].idxmin()
                # best_compression_tool.append(data['compression'].iloc[idx])

            # plot methods vs sizes
            fig = plt.figure()
            plt.title(f"{sequence} - {counts_name} - k={kmer_size}")
            plt.xlabel("methods")
            plt.ylabel("sizes [bytes]")
            plt.bar(methods, no_compression_size)
            plt.bar(methods, best_compression_size)
            addlabels(methods, no_compression_size, best_compression_size)
            plt.savefig(f"{FIGURES_PATH}/{sequence}.{counts_name}.k{kmer_size}.png")
            plt.clf()

            print(f"\t\tuncompressed: {no_compression_size}")
            print(f"\t\tcompressed: {best_compression_size}")


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
    plot_sequences_vs_methods(sequences, kmer_sizes, data, False)
    # counts
    # plot_sequences_vs_methods(sequences, kmer_sizes, data, True)


if __name__ == "__main__":
    main()
