import os
import argparse
import methods

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

parser.add_argument("-sc", "--snakemake_cores", type=str,
                    help="Number of cores for Snakemake", default=6)

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

args = parser.parse_args()

rgi = args.annotations
reads_path = args.reads
envs_path = args.envs
metadata = args.metadata
sc = args.snakemake_cores

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
rp_order = methods.nsort(reads_path, False)
rp_path_order = methods.nsort(reads_path, True)
rpnum_rp = {}

#create kallisto and shortbred folders
# os.mkdir(f"{outdir}/kallisto")
# os.mkdir(f"{outdir}/shortbred")

#create the FASTA and .faa files
uid_tracker, protein_tracker, fasta_path, faa_path, head = methods.fasta(f"{outdir}/all_rgi", f"{outdir}/fasta")

#create config file for abundance run
d = {"output": outdir, "reads": reads_path, "sample": "N/A",
    "conda_path": conda_profile, "envs_path": envs_path,
    "illuminaclip": "N/A", "all_args": fasta_path,
    "muscle": f"{dep}/muscle3.8.31_i86linux64", 
    "usearch": f"{dep}/usearch11.0.667_i86linux32", 
    "CARD_markers": f"{dep}/ShortBRED_CARD_2017_markers.faa", "rule_all": "abundance"}

#create 2nd config file
config_path = methods.config(outdir, d)

#second call of snakemake
os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir} all --configfile {config_path}")

#gather all kallisto and shortbred output files
os.system(f"bash {home_dir}/scripts/gather_abundance.sh {rp_total-1} {outdir}")

#gather kallisto and shortbred counts
first_col = methods.counts(uid_tracker, f"{outdir}/all_kallisto", f"{outdir}/all_shortbred", f"{outdir}/matrices")

#create observation table
methods.observations(uid_tracker, head, first_col, metadata, f"{outdir}/matrices")