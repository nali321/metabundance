#!/bin/bash

mkdir $2/all_snp

for ((i=1; i<=$1; i++)); do
    cp $2/spraynpray/read_pair_"$i"/read_pair_"$i".csv $2/all_snp/"$i"_snp.csv
done