import os
import argparse
import methods

#pipeline steps:
#PART 1:
#1. trim read pairs
#2. assemble metagenome
#3. annotate metagenome with RGI - DO STEPS 1-3 SEPARATELY

#PART 2:
#4. Data processing steps: - START HERE FOR TESTING PURPOSES!!!
    #a. Create a FASTA file of all DNA sequences from RGI
    #b. Create a protein FASTA file for each sample with only ARGs from it
#5. Run Kallisto by mapping read pairs back to ARG FASTA file
#6. Create input files for Metagenomseq:
    #a. Create ARG - read pair counts matrix
    #b. Create taxa table where OTU 1-Max ARG Number is matched with it's subsequent RGI data row
    #c. Supply phenotype file from user
#7. Create Metagenomseq object in R and obtain normalized counts matrix
#8. Use shortbred identify and quantify with individual ARG .faa files
#9. Create files for phyloseq to use

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--annotations", type=str,
                    help="Filepath to folder with ARG annotations", required=True)

parser.add_argument("-r", "--read_pairs", type=str,
                    help="Filepath to folder of read pairs", required=True)

parser.add_argument("-e", "--envs", type=str,
                    help="Filepath to environments folder", required=True)

parser.add_argument("-m", "--metadata", type=str,
                    help="Metadata file. First column must be ids of each sample. \
                        Rest of columns can be phenotypic information about each sample", required=True)

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

# parser.add_argument("-sc", "--snakemake_cores", type=str,
#                     help="Number of cores for Snakemake", default=4)

args = parser.parse_args()

rgi = args.annotations
read_pairs = args.read_pairs
envs_path = args.envs
metadata = args.metadata

#create output directory
outdir = args.outdir
try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#get the snakefile and dependencies path
script_dir = os.path.dirname(os.path.abspath(__file__))
snake_dir = os.path.join(script_dir, "workflow/Snakefile")
dep = os.path.join(envs_path, "dependencies").replace("\\", "/")

#get conda profile path
conda_path_file = os.path.join(outdir, "CONDA_PATH.txt").replace("\\", "/")
os.system(f"conda info --root > {conda_path_file}")
with open (conda_path_file, 'r') as file:
    for line in file:
        conda_profile = os.path.join(line.strip('\n'), "etc/profile.d/conda.sh").replace("\\", "/")
        break

#create the FASTA and .faa files
fasta_dir = os.path.join(rgi, f"{outdir}/fasta")
uid_tracker, protein_tracker, fasta_path, faa_path, head = methods.fasta(rgi, fasta_dir)

# #run the kallisto index rule
# os.system(f"snakemake --directory {outdir} --snakefile {snake_dir} kallisto_index -c6 --configfile {config_path}")

#order read pairs folder
rp_order = methods.nsort(read_pairs, False)
rp_path_order = methods.nsort(read_pairs, True)
rpnum_rp = {}

#create kallisto and shortbred folders
# os.mkdir(f"{outdir}/kallisto")
# os.mkdir(f"{outdir}/shortbred")

#go through each read pair and perform kallisto + shortbred on it
rpnum = 1
for i in range(0, len(rp_order), 2):
    rp1 = rp_path_order[i]
    rp2 = rp_path_order[i+1]
    # rp2 = rp_path_order[i+1] if i+1 < len(rp_order) else None
    # rpnum = rp_order[i].split("_")[0]
    
    #create config dictionary
    d = {"output": outdir, "forward": rp1, "reverse": rp2, "conda_path": conda_profile,
        "envs_path": envs_path, "illuminaclip": "N/A", "all_args": f"{fasta_path}",
        "id": rpnum, "muscle": f"{dep}/muscle3.8.31_i86linux64", "usearch": f"{dep}/usearch11.0.667_i86linux32", 
        "CARD_markers": f"{dep}/ShortBRED_CARD_2017_markers.faa", "faa_file": f"{faa_path}/{rpnum}.faa",
        "targets": ["kallisto", "shortbred_quantify"], "target_output": {"kallisto": f"{outdir}/kallisto/read_pair_{rpnum}/abundance.tsv", 
                                                                         "shortbred_quantify": f"{outdir}/shortbred/shortbred_quantify/read_pair_{rpnum}/results.tsv"}}

    #create config file
    config_path = methods.config(outdir, d)

    #call snakemake with 4 cores and call two rules, expecting 4 rules to run at once in parallel
    # os.system(f"snakemake -j 4 --directory {outdir} --snakefile {snake_dir} kallisto shortbred_quantify -c6 --configfile {config_path}")
    os.system(f"snakemake -j 4 --directory {outdir} --snakefile {snake_dir} -c6 --configfile {config_path}")

    #match rpnum to rp filename
    rpnum_rp[rpnum] = rp_order[i]
    rpnum +=1

#create folder collection and otu folder
os.mkdir(f"{outdir}/kallisto/kallisto_files")
os.mkdir(f"{outdir}/shortbred/shortbred_files")
os.mkdir(f"{outdir}/otus")

#collect all kallisto and shortbred files
os.system(f"bash {script_dir}/gather.sh {rpnum-1} {outdir}")

#gather kallisto and shortbred counts
first_col = methods.counts(uid_tracker, f"{outdir}/kallisto/kallisto_files", f"{outdir}/shortbred/shortbred_files", f"{outdir}/otus")

#create otu tables
methods.otu(uid_tracker, head, first_col, metadata, f"{outdir}/otus")