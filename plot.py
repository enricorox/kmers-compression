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
def addlabels2(x, y1, y2):
    for i in range(len(x)):
        plt.text(i, y1[i], f"{y1[i]}", ha='center', va='bottom', color='green')
        plt.text(i, y2[i], f"{y2[i]}", ha='center', va='top', color='green')


def addlabels1(x, y):
    for i in range(len(x)):
        plt.text(i, y[i], f"{y[i]:.2}", ha='center', va='bottom', color='green')


def plot(sequences: list, kmer_sizes: list, data: pd.DataFrame, counts: bool):
    if counts:
        counts_name = "counts"
        methods = [e for e in list(METHODS) if e not in NO_COUNTS_METHOD]
    else:
        counts_name = "no-counts"
        methods = list(METHODS)
    print(counts_name)
    all_sizes = {}
    for kmer_size in kmer_sizes:
        all_sequences = {}
        print(f"\tkmer-size: {kmer_size}")
        avg_ratios = [0 for i in range(len(methods))]
        for sequence in sequences:
            print(f"\t\tsequence: {sequence}")
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
            for i, method in enumerate(methods):
                # sum all sizes of uncompressed files for this method
                size = seq_df.query(f'method == "{method}" and compression == "none"')['size'].sum()
                no_compression_size.append(size)

                # find the best compression size
                comp_size = 0
                # consider all file types
                for file_type in METHODS.get(method):
                    q = seq_df.query(f'method == "{method}" and compression != "none" and file_type == "{file_type}"')
                    comp_size += q['size'].min()
                    # only one file if no counts
                    if not counts:
                        break
                best_compression_size.append(comp_size)

                # find the best compression tool
                # idx = q['size'].idxmin()
                # best_compression_tool.append(data['compression'].iloc[idx])

                avg_ratios[i] += (size / comp_size) / len(methods)
            print(f"\t\t\tuncompressed: {no_compression_size}")
            print(f"\t\t\tcompressed: {best_compression_size}")

            # plot methods vs sizes
            plt.title(f"{sequence} - {counts_name} - k={kmer_size}")
            plt.xlabel("methods")
            plt.ylabel("sizes [bytes]")
            plt.bar(methods, no_compression_size)
            plt.bar(methods, best_compression_size)
            addlabels2(methods, no_compression_size, best_compression_size)
            plt.savefig(f"{FIGURES_PATH}/{sequence}.{counts_name}.k{kmer_size}.png")
            plt.clf()

            # save sizes
            all_sequences.update({sequence: [no_compression_size, best_compression_size]})

        # plot average rations
        plt.title(f"average - {counts_name} - k={kmer_size}")
        plt.xlabel("methods")
        plt.ylabel("size ratios")
        plt.bar(methods, avg_ratios)
        addlabels1(methods, avg_ratios)
        plt.savefig(f"{FIGURES_PATH}/average.{counts_name}.k{kmer_size}.png")
        plt.clf()
        all_sizes.update({kmer_size: all_sequences})
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
    print(data.at[0, 'size'])
    exit(1)
    # read sequence file
    sequences = read_from_file(SEQUENCES_FILE)

    # read kmer-sizes file
    kmer_sizes_file = "kmer-sizes.txt"
    kmer_sizes = read_from_file(kmer_sizes_file)

    # no-counts
    plot(sequences, kmer_sizes, data, False)
    # counts
    plot(sequences, kmer_sizes, data, True)


if __name__ == "__main__":
    main()
