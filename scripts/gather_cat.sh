#!/bin/bash

mkdir $2/all_cat

for ((i=1; i<=$1; i++)); do
    cp $2/cat/read_pair_"$i"/CAT_"$i".contig2classification.txt $2/all_cat/"$i"_cat.csv
done