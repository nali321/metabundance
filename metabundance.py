import os
import argparse
import methods

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reads", type=str,
                    help="Filepath to reads", required=True)

parser.add_argument("-e", "--envs", type=str,
                    help="Filepath to environments folder", required=True)

parser.add_argument("-sc", "--snakemake_cores", type=str,
                    help="Number of cores for Snakemake", default="6")

parser.add_argument("-i", "--illuminaclip", type=str,
                    help="Choose the Illuminaclip that will be used with Trimmomatic", required=True)

parser.add_argument("-m", "--metadata", type=str,
                    help="Metadata file. First column must be ids of each sample. \
                        Rest of columns can be phenotypic information about each sample", required=True)

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

args = parser.parse_args()

reads_path = args.reads
envs_path = args.envs
sc = args.snakemake_cores
illuminaclip = args.illuminaclip
metadata = args.metadata

#create output directory
outdir = args.outdir
try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#get the snakefile, script, and dependencies path
home_dir = os.path.dirname(os.path.abspath(__file__))
snake_dir = os.path.join(home_dir, "workflow/Snakefile")
scripts_dir = os.path.join(home_dir, "workflow/scripts")
dep = os.path.join(envs_path, "dependencies").replace("\\", "/")

#get conda profile path
conda_path_file = os.path.join(outdir, "CONDA_PATH.txt").replace("\\", "/")
os.system(f"conda info --root > {conda_path_file}")
with open (conda_path_file, 'r') as file:
    for line in file:
        conda_profile = os.path.join(line.strip('\n'), "etc/profile.d/conda.sh").replace("\\", "/")
        break

#prepare reads for snakemake use
rp_order = methods.nsort(reads_path, True)

#check if folder is good
methods.check_reads(len(rp_order))

#get total number of read pairs
rp_total = int(len(rp_order)/2)

#get sample ids
sample_ids = [i for i in range(1, rp_total+1)]

#create config file for rgi run
d = {"output": outdir, "reads": reads_path, "sample": sample_ids,
    "conda_path": conda_profile, "envs_path": envs_path,
    "illuminaclip": illuminaclip, "fasta": "N/A",
    "muscle": "N/A", "usearch": "N/A", "CARD_markers": "N/A", "rule_all": "annotations"}

#create config file
config_path = methods.config(outdir, d)

#call snakemake
os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir} genomad rgi integron_finder --configfile {config_path}")

# #collect all rgi, genomad, integron output files
# os.system(f"bash {home_dir}/scripts/gather_annotations.sh {rp_total-1} {outdir}")

# #create the FASTA and .faa files
# uid_tracker, protein_tracker, fasta_path, faa_path, head = methods.fasta(f"{outdir}/all_rgi", f"{outdir}/fasta")

# #create config file for abundance run
# d = {"output": outdir, "reads": reads_path, "sample": sample_ids,
#     "conda_path": conda_profile, "envs_path": envs_path,
#     "illuminaclip": illuminaclip, "all_args": fasta_path,
#     "muscle": f"{dep}/muscle3.8.31_i86linux64", 
#     "usearch": f"{dep}/usearch11.0.667_i86linux32", 
#     "CARD_markers": f"{dep}/ShortBRED_CARD_2017_markers.faa", "rule_all": "abundance"}

# #create 2nd config file
# config_path = methods.config(outdir, d)

# #second call of snakemake
# os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir} all --configfile {config_path}")

# #gather all kallisto and shortbred output files
# os.system(f"bash {home_dir}/scripts/gather_abundance.sh {rp_total-1} {outdir}")

# #gather kallisto and shortbred counts
# first_col = methods.counts(uid_tracker, f"{outdir}/all_kallisto", f"{outdir}/all_shortbred", f"{outdir}/matrices")

# #create observation table
# methods.observations(uid_tracker, head, first_col, metadata, f"{outdir}/matrices")