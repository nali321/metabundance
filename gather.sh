#!/bin/bash

for ((i=1; i<=$1; i++)); do
    cp "$2"/kallisto/read_pair_"$i"/abundance.tsv "$2"/kallisto/kallisto_files/"$i"_abundance.tsv
    cp "$2"/shortbred/shortbred_quantify/read_pair_"$i"/results.tsv "$2"/shortbred/shortbred_files/"$i"_results.tsv
done