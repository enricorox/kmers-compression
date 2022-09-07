#!/bin/env python3
import sys

import pandas as pd
import matplotlib.pyplot as plt

defalut_results = "results/results_2022-09-07_22.44.50.csv"
sequences_file = "sequences-test.txt"
methods = ["bcalm", "ust", "prophasm", "metagraph", "assembly"]
figures_path = "figures"


def read_from_file(myfile):
    seq = []
    with open(myfile, "r") as file:
        for line in file:
            sequence = line.strip()
            if sequence == "" or sequence[0] == "#":
                continue
            seq.append(sequence)
    return seq


def main():
    print("*** Executing plot utility...")
    # parse parameter or set default
    if len(sys.argv) == 2:
        infile = sys.argv[1]
    else:
        infile = defalut_results

    # read csv file
    # headers: sequence,method,counts,kmer_size,file_type,compression,size
    data = pd.read_csv(infile)
    print(data)
    # read sequence file
    sequences = read_from_file(sequences_file)

    # read kmer-sizes file
    kmer_sizes_file = "kmer-sizes.txt"
    kmer_sizes = read_from_file(kmer_sizes_file)

    # no-counts
    for kmer_size in kmer_sizes:
        for sequence in sequences:
            # select sequence and kmer-size
            seq_df = data.query(
                f'sequence == "{sequence}" and counts == "no-counts"'
            ).query(f'kmer_size == {kmer_size}')
            print(f"{sequence},{kmer_size}")
            print(seq_df)
            exit()
            # for each method find the size
            no_compression = []
            for method in methods:
                size = seq_df.query(f'method == "{method}" and compression == "none"')['size'].iloc[0]
                no_compression.append(size)

            # for each method find the compressed size and best compression tool
            best_compression_size = []
            best_compression_tool = []
            for method in methods:
                size = seq_df.query(f'method == "{method}" and compression != "none"')['size'].min()
                idx = seq_df.query(f'method == "{method}" and compression != "none"')['size'].idxmin()
                best_compression_size.append(size)
                best_compression_tool.append(data['compression'].iloc[idx])

            # plot
            fig = plt.figure()
            plt.title(f"{sequence} - k={kmer_size}")
            plt.bar(methods, no_compression)
            plt.bar(methods, best_compression_size)
            plt.savefig(f"{figures_path}/{sequence}.k{kmer_size}.png")



if __name__ == "__main__":
    main()
