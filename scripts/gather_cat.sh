#!/bin/bash

mkdir $2/all_cat

for ((i=1; i<=$1; i++)); do
        # Copy CAT file
    if [[ -f "$2/cat/read_pair_$i/CAT_$i.contigs2classification.txt" ]]; then
        cp "$2/cat/read_pair_$i/CAT_$i.contigs2classification.txt" "$2/all_cat/${i}_cat.txt"
    else
        echo "CAT file for read_pair_$i does not exist. Skipping..."
    fi
done