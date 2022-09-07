#!/bin/env python3

import pandas as pd

# read csv file
infile = "results/results_2022-09-07_16.30.08.csv"
data = pd.read_csv(infile)

# read sequence file
sequences_file = "sequences-test.txt"
sequences = []
with open(sequences_file, "r") as file:
    next(file)
    for line in file:
        sequence = line.strip(" \n")
        if sequence == "" or sequence[0] == "#":
            continue
        sequences.append(sequence)

print(sequences)
f = data.query('method == "metagraph"')

#print(f)
