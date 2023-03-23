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

parser.add_argument("-1", "--read_1", type=str,
                    help="Filepath to forward read", required=True)

parser.add_argument("-2", "--read_2", type=str,
                    help="Filepath to reverse read", required=True)

parser.add_argument("-e", "--envs", type=str,
                    help="Filepath to environments folder", required=True) 

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

# parser.add_argument("-sc", "--snakemake_cores", type=str,
#                     help="Number of cores for Snakemake", default=4)

parser.add_argument("-i", "--illuminaclip", type=str,
                    help="Choose the Illuminaclip that will be used with Trimmomatic", required=True)

args = parser.parse_args()

r1 = args.read_1
r2 = args.read_2
envs_path = args.envs
illuminaclip = args.illuminaclip

#create output directory
outdir = args.outdir
try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#get the snakefile path
script_dir = os.path.dirname(os.path.abspath(__file__))
snake_dir = os.path.join(script_dir, "workflow/Snakefile")

#get conda profile path
conda_path_file = os.path.join(outdir, "CONDA_PATH.txt").replace("\\", "/")
os.system(f"conda info --root > {conda_path_file}")
with open (conda_path_file, 'r') as file:
    for line in file:
        conda_profile = os.path.join(line.strip('\n'), "etc/profile.d/conda.sh").replace("\\", "/")
        break

#create config file
d = {"output": outdir, "forward": r1, "reverse": r2, "conda_path": conda_profile,
     "envs_path": envs_path, "illuminaclip": illuminaclip, "all_args": "N/A",
     "id": "N/A", "muscle": "N/A", "usearch": "N/A", "CARD_markers": "N/A", "faa_file": "N/A"}

#create config file
config_path = methods.config(outdir, d)

os.system(f"snakemake --directory {outdir} --snakefile {snake_dir} rgi -c6 --configfile {config_path}")
# os.system(f"snakemake --snakefile /mmfs1/home/4565alin/READS2ARGS/workflow/Snakefile -c6 --use-conda --conda-create-envs-only")

#potentially add onto a file in the main output folder to create a read pair - assembly dictionary
#or even add onto it with the "final.txt" file from rgi