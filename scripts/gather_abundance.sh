#!/bin/bash

mkdir $2/all_kallisto
mkdir $2/all_shortbred

for ((i=1; i<=$1; i++)); do
    cp $2/kallisto/read_pair_"$i"/abundance.tsv $2/all_kallisto/"$i"_abundance.tsv
    cp $2/shortbred_quantify/read_pair_"$i"/results.tsv $2/all_shortbred/"$i"_results.tsv
done