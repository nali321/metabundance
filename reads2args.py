import os
import argparse
import methods


parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reads", type=str,
                    help="Filepath to reads", required=True)

parser.add_argument("-e", "--envs", type=str,
                    help="Filepath to environments folder", required=True) 

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where output will go", required=True)

parser.add_argument("-sc", "--snakemake_cores", type=str,
                    help="Number of cores for Snakemake", default=6)

parser.add_argument("-i", "--illuminaclip", type=str,
                    help="Choose the Illuminaclip that will be used with Trimmomatic", required=True)

args = parser.parse_args()

reads_path = args.reads
envs_path = args.envs
sc = args.sc
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

os.system(f"snakemake --directory {outdir} --snakefile {snake_dir} rgi -c{sc} --configfile {config_path}")
# os.system(f"snakemake --snakefile /mmfs1/home/4565alin/READS2ARGS/workflow/Snakefile -c6 --use-conda --conda-create-envs-only")

#potentially add onto a file in the main output folder to create a read pair - assembly dictionary
#or even add onto it with the "final.txt" file from rgi