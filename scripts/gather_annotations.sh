#!/bin/bash

mkdir $2/all_rgi
mkdir $2/all_plasmid
mkdir $2/all_virus
mkdir $2/all_integron
mdkir $2/all_ice
mkdir $2/all_tn_is

for ((i=1; i<=$1; i++)); do
    # Copy RGI file
    if [[ -f "$2/rgi/read_pair_$i/final.txt" ]]; then
        cp "$2/rgi/read_pair_$i/final.txt" "$2/all_rgi/${i}_rgi.txt"
    else
        echo "RGI file for read_pair_$i does not exist. Skipping..."
    fi

    # Copy plasmid summary
    if [[ -f "$2/genomad/read_pair_$i/contigs_summary/contigs_plasmid_summary.tsv" ]]; then
        cp "$2/genomad/read_pair_$i/contigs_summary/contigs_plasmid_summary.tsv" "$2/all_plasmid/${i}_plasmid.tsv"
    else
        echo "Plasmid summary for read_pair_$i does not exist. Skipping..."
    fi

    # Copy virus summary
    if [[ -f "$2/genomad/read_pair_$i/contigs_summary/contigs_virus_summary.tsv" ]]; then
        cp "$2/genomad/read_pair_$i/contigs_summary/contigs_virus_summary.tsv" "$2/all_virus/${i}_virus.tsv"
    else
        echo "Virus summary for read_pair_$i does not exist. Skipping..."
    fi

    # Copy integrons file
    if [[ -f "$2/integron_finder/read_pair_$i/Results_Integron_Finder_contigs/contigs.integrons" ]]; then
        cp "$2/integron_finder/read_pair_$i/Results_Integron_Finder_contigs/contigs.integrons" "$2/all_integron/${i}_integrons.txt"
    else
        echo "Integron file for read_pair_$i does not exist. Skipping..."
    fi

    # Copy ICE file
    if [[ -f "$2/ice/read_pair_$i/ice_$i.txt" ]]; then
        cp "$2/ice/read_pair_$i/ice_$i.txt" "$2/all_ice/${i}_ice.txt"
    else
        echo "ICE file for read_pair_$i does not exist. Skipping..."
    fi

    # Copy TN/IS file
    if [[ -f "$2/tn_is/read_pair_$i/tn_is_$i.txt" ]]; then
        cp "$2/tn_is/read_pair_$i/tn_is_$i.txt" "$2/all_tn_is/${i}_tn_is.txt"
    else
        echo "TN/IS file for read_pair_$i does not exist. Skipping..."
    fi
done