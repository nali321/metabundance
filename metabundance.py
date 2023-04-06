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
snake_dir = os.path.join(home_dir, "workflow")
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
    "illuminaclip": illuminaclip, "fasta": "N/A", "protein": "N/A", "snp": "N/A",
    "muscle": "N/A", "usearch": "N/A", "CARD_markers": "N/A", "rule_all": "annotations"}

#create config file
config_path = methods.config(d, "config1", outdir)

#call annotations snakefile
os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir}/Snakefile all --configfile {config_path}")

#collect all rgi, genomad, integron output files
os.system(f"bash {home_dir}/scripts/gather_annotations.sh {rp_total} {outdir}")

#create the FASTA and .faa files
#if you're re-running the snakefile and the annotations.FASTA or .faa files get re-made due to
#this being re-ran, then it'll re-run the kallisto and shortbred rules because the input files
#were updated since the previous run
# uid_tracker, protein_tracker, fasta_path, faa_path, head = methods.fasta(f"{outdir}/all_rgi", f"{outdir}/fasta")
uid_tracker, protein_tracker, head = methods.fasta(f"{outdir}/all_rgi", f"{outdir}/fasta")

#create FASTA file
fasta_path = f"{outdir}/fasta"
try:
    os.mkdir(outdir)
#throw an error if the folder doesn't exist
except OSError as error:
    print(error)
#if the folder was just created, create the input files. this way you can re-run snakefile and have the
#kallisto and shortbred steps not be re-ran
else:
    with open (os.path.join(fasta_path, "annotations.FASTA").replace("\\", "/"), 'w+') as f:
        for x in uid_tracker:
            f.write(f">{x}|{uid_tracker[x][0]}\n{uid_tracker[x][1][17]}\n")
            if uid_tracker[x][0] not in protein_tracker:
                protein_tracker[uid_tracker[x][0]] = [[x, uid_tracker[x][1], uid_tracker[x][2]]]
            else:
                protein_tracker[uid_tracker[x][0]].append([x, uid_tracker[x][1], uid_tracker[x][2]])

#create .faa files
faa_path = os.path.join(fasta_path, "protein_files")
try:
    os.mkdir(os.path.join(faa_path, "protein_files").replace("\\", "/"))
except OSError as error:
    print(error)
else:
    for x in protein_tracker:
        with open (os.path.join(faa_path, f"protein_files/{x}.faa").replace("\\", "/"), 'w+') as f:
            for y in protein_tracker[x]:
                f.write(f">{y[0]}|{x}\n{y[2]}\n")

#create config file for abundance run
d = {"output": outdir, "reads": reads_path, "sample": sample_ids,
    "conda_path": conda_profile, "envs_path": envs_path,
    "illuminaclip": illuminaclip, "fasta": f"{fasta_path}/annotations.FASTA",
    "protein": faa_path, "snp": "N/A", "muscle": f"{dep}/muscle3.8.31_i86linux64",
    "usearch": f"{dep}/usearch11.0.667_i86linux32", 
    "CARD_markers": f"{dep}/ShortBRED_CARD_2017_markers.faa", "rule_all": "abundance"}

#create 2nd config file
config_path2 = methods.config(d, "config2", outdir)

#second call of snakemake
os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir}/Snakefile all --configfile {config_path2}")

#gather all kallisto and shortbred output files
os.system(f"bash {home_dir}/scripts/gather_abundance.sh {rp_total} {outdir}")

#gather kallisto and shortbred counts
first_col = methods.counts(uid_tracker, f"{outdir}/all_kallisto", f"{outdir}/all_shortbred", f"{outdir}/matrices")

#create observation table
methods.observations(uid_tracker, head, first_col, metadata, f"{outdir}/matrices")

#attach genomad, integron, and spraynpray (if needed)