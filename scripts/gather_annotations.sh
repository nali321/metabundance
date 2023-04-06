#!/bin/bash

mkdir $2/all_rgi
mkdir $2/all_plasmid
mkdir $2/all_virus
mkdir $2/all_integron

for ((i=1; i<=$1; i++)); do
    cp $2/rgi/read_pair_"$i"/final.txt $2/all_rgi/"$i"_rgi.txt
    cp $2/genomad/read_pair_"$i"/contigs_summary/contigs_plasmid_summary.tsv $2/all_plasmid/"$i"_plasmid.tsv
    cp $2/genomad/read_pair_"$i"/contigs_summary/contigs_virus_summary.tsv $2/all_virus/"$i"_virus.tsv
    cp $2/integron_finder/read_pair_"$i"/Results_Integron_Finder_contigs/contigs.integrons" $2/all_integron/"$i".integrons
done