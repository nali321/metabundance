import os
import argparse
import functions

parser = argparse.ArgumentParser()

#required flags

parser.add_argument("-r", "--reads", type=str,
                    help="Filepath to reads", required=True)

parser.add_argument("-e", "--envs", type=str,
                    help="Filepath to environments folder", required=True)

parser.add_argument("-i", "--illuminaclip", type=str,
                    help="Choose the Illuminaclip that will be used with Trimmomatic", required=True)

parser.add_argument("-m", "--metadata", type=str,
                    help="Metadata file. First column must be ids of each sample. \
                        Rest of columns can be phenotypic information about each sample", required=True)

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

#optional flags

parser.add_argument("-sc", "--snakemake_cores", type=str,
                    help="Number of cores for Snakemake to use. Default is 6", default="6")

#these three flags must be submitted together

parser.add_argument("-cp", "--cat_path", type=str,
                    help="Include the filepath to the CAT conda environment if you want taxonomic identification", default="N/A")

parser.add_argument("-cpk", "--cat_pack", type=str,
                    help="Include the filepath to the CAT pack installation if you want taxonomic identification", default="N/A")

parser.add_argument("-cd", "--cat_db", type=str,
                    help="Include the filepath to the CAT database folder if you want taxonomic identification (not the db itself but the folder that the db is in)", default="N/A")

args = parser.parse_args()

reads_path = args.reads
envs_path = args.envs
sc = args.snakemake_cores
illuminaclip = args.illuminaclip
cat_path = args.cat_path
cat_db = args.cat_db
cat_pack = args.cat_pack
metadata = args.metadata

#check if taxonomy needs to be ran
if cat_path and cat_pack and cat_db == "N/A":
    rule_type = "annotations"
elif cat_path and cat_pack and cat_db != "N/A":
    rule_type = "taxonomy"
#throw errors if only one or two of the paths is provided
na_count = sum(var == "N/A" for var in [cat_path, cat_pack, cat_db])
if na_count != 0 and na_count < 3:
    raise Exception('All three filepaths to cat_path, cat_pack, and cat_db need to be provided')

# elif cat_path == "N/A" and cat_db != "N/A":
#     raise Exception('Both filepaths to cat_path and cat_db need to be provided')
# elif cat_path != "N/A" and cat_db == "N/A":
#     raise Exception('Both filepaths to cat_path and cat_db need to be provided')

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
rp_order = functions.nsort(reads_path, True)

#check if folder is good
functions.check_reads(len(rp_order))

#get total number of read pairs
rp_total = int(len(rp_order)/2)

#get sample ids
sample_ids = [i for i in range(1, rp_total+1)]

#obtain conda profile from cat path in case its from a different conda installation
parts = cat_path.split("miniconda3", 1)
extracted_path = parts[0] + "miniconda3"
cat_profile = f"{extracted_path}/etc/profile.d/conda.sh"

#create config file for rgi run
d = {"output": outdir, "reads": reads_path, "sample": sample_ids, "conda_path": conda_profile, "envs_path": envs_path,
    "illuminaclip": illuminaclip, "fasta": "N/A", "protein": "N/A", "cat_path": cat_path, "cat_db": cat_db, "cat_profile": cat_profile, "cat_pack": cat_pack,
    "ice_db": f"{dep}/tncentral_isfinder/tncentral_isfinder.fa", "tn_is_db": f"{dep}/ICEberg/ICE_seq_all.fas", "muscle": "N/A", "usearch": "N/A", "CARD_markers": "N/A", "rule_all": "DEBUG"}

#create config file
config_path = functions.config(d, "config1", outdir)

#call annotations snakefile
os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir}/Snakefile all --configfile {config_path}")

# #pipeline ending early here

# #collect all annotationj output files
# os.system(f"bash {home_dir}/scripts/gather_annotations.sh {rp_total} {outdir}")

# #create the FASTA and .faa files
# #if you're re-running the snakefile and the annotations.FASTA or .faa files get re-made due to
# #this being re-ran, then it'll re-run the kallisto and shortbred rules because the input files
# #were updated since the previous run
# # uid_tracker, protein_tracker, fasta_path, faa_path, head = methods.fasta(f"{outdir}/all_rgi", f"{outdir}/fasta")
# uid_tracker, protein_tracker, head = functions.fasta(f"{outdir}/all_rgi")

# #create FASTA file
# fasta_path = f"{outdir}/fasta"
# try:
#     os.mkdir(fasta_path)
# #throw an error if the folder already exists
# except OSError as error:
#     print(error)
# #if the folder was just created, create the input files. this way you can re-run snakefile and have the
# #kallisto and shortbred steps not be re-ran
# else:
#     with open (os.path.join(fasta_path, "annotations.FASTA").replace("\\", "/"), 'w+') as f:
#         for x in uid_tracker:
#             f.write(f">{x}|{uid_tracker[x][0]}\n{uid_tracker[x][1][17]}\n")
#             if uid_tracker[x][0] not in protein_tracker:
#                 protein_tracker[uid_tracker[x][0]] = [[x, uid_tracker[x][1], uid_tracker[x][2]]]
#             else:
#                 protein_tracker[uid_tracker[x][0]].append([x, uid_tracker[x][1], uid_tracker[x][2]])

# #create .faa files
# faa_path = os.path.join(fasta_path, "protein_files")
# try:
#     os.mkdir(faa_path)
# except OSError as error:
#     print(error)
# else:
#     for x in protein_tracker:
#         with open (os.path.join(faa_path, f"{x}.faa").replace("\\", "/"), 'w+') as f:
#             for y in protein_tracker[x]:
#                 f.write(f">{y[0]}|{x}\n{y[2]}\n")

# #create config file for abundance run
# d = {"output": outdir, "reads": reads_path, "sample": sample_ids,
#     "conda_path": conda_profile, "envs_path": envs_path,
#     "illuminaclip": illuminaclip, "fasta": f"{fasta_path}/annotations.FASTA",
#     "protein": faa_path, "cat_path": "N/A", "cat_db": "N/A", "cat_profile": "N/A", "cat_pack": "N/A", "ice_db": "N/A", "tn_is_db": "N/A",
#     "muscle": f"{dep}/muscle3.8.31_i86linux64", "usearch": f"{dep}/usearch11.0.667_i86linux32", 
#     "CARD_markers": f"{dep}/ShortBRED_CARD_2017_markers.faa", "rule_all": "abundance"}

# #create 2nd config file
# config_path2 = functions.config(d, "config2", outdir)

# #second call of snakemake
# os.system(f"snakemake --cores {sc} --directory {outdir} --snakefile {snake_dir}/Snakefile all --configfile {config_path2}")

# #gather all kallisto and shortbred output files
# os.system(f"bash {home_dir}/scripts/gather_abundance.sh {rp_total} {outdir}")

# #gather kallisto and shortbred counts
# first_col = functions.count_matrices(uid_tracker, f"{outdir}/all_kallisto", f"{outdir}/all_shortbred", f"{outdir}/matrices")

# #create observation table
# functions.observations(uid_tracker, head, first_col, metadata, f"{outdir}/matrices")

# #attach MGE and CAT (if needed) data after the observations matrix has been made
# functions.genomad(f"{outdir}/all_plasmid", f"{outdir}/all_virus", f"{outdir}/matrices")
# functions.blast_mges(f"{outdir}/all_ice", f"{outdir}/matrices")
# functions.blast_mges(f"{outdir}/all_tn_is", f"{outdir}/matrices")
# functions.integron(f"{outdir}/all_integron", f"{outdir}/matrices")

# #gather snp files and append data to observation matrix if necessary
# if rule_type == "taxonomy":
#     os.system(f"bash {home_dir}/scripts/gather_cat.sh {rp_total} {outdir}")
#     functions.cat(f"{outdir}/all_cat", f"{outdir}/matrices")