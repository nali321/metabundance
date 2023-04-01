import os
import argparse

parser = argparse.ArgumentParser()

# parser.add_argument("-c", "--conda_profile", type=str,
#                     help="Filepath to conda profile", required=True)

parser.add_argument("-o", "--outdir", type=str,
                    help="Directory where environments will go", required=True)

args = parser.parse_args()

#create output directory
outdir = args.outdir
try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#get conda profile path
conda_path_file = os.path.join(outdir, "CONDA_PATH.txt").replace("\\", "/")
os.system(f"conda info --root > {conda_path_file}")
with open (conda_path_file, 'r') as file:
    for line in file:
        conda_profile = os.path.join(line.strip('\n'), "etc/profile.d/conda.sh").replace("\\", "/")
        break

#get the script path
home_dir = os.path.dirname(os.path.abspath(__file__))
envs_dir = os.path.join(home_dir, "workflow/envs")
usearch_dir = os.path.join(home_dir, "workflow/dependencies/usearch11.0.667_i86linux32.gz")
muscle_dir = os.path.join(home_dir, "workflow/dependencies/muscle3.8.31_i86linux64.tar.gz")
cardmarkers_dir = os.path.join(home_dir, "workflow/dependencies/ShortBRED_CARD_2017_markers.faa.gz")

#create switches
blast_switch = os.path.join(home_dir, "scripts/blast_switch.sh").replace("\\", "/")
genomad_switch = os.path.join(home_dir, "scripts/genomad_switch.sh").replace("\\", "/")

#iterate over predownloaded envs folder, and download them from the .yamls
for filename in os.listdir(envs_dir):
    #create environment from .yaml file
    # if filename != "blast.yaml":
    if filename == "genomad.yaml":
        z = os.path.join(envs_dir, filename).replace("\\", "/")
        env = os.path.join(outdir, filename.split(".")[0])
        os.mkdir(env)
        os.system(f"mamba env create -f {z} -p {env}")
    
    # #blast cannot be created from a .yaml file, download it as a separate env
    # if filename == "blast.yaml":
    #     env = os.path.join(outdir, "blast_env")
    #     os.system(f"bash {blast_switch} {conda_profile} {env}")
    
    # if filename == "shortbred.yaml":
    #     #create dependencies folder and cd into it
    #     dep = os.path.join(outdir, "dependencies").replace("\\", "/")
    #     os.mkdir(dep)
    #     os.chdir(dep)

    #     #install muscle
    #     os.system(f"cp {muscle_dir} {dep}")
    #     os.system(f"tar -zxf muscle3.8.31_i86linux64.tar.gz")

    #     #install usearch
    #     os.system(f"cp {usearch_dir} {dep}")
    #     os.system(f"gzip -d usearch11.0.667_i86linux32.gz")
    #     os.system(f"chmod +x usearch11.0.667_i86linux32")

    #     #unzip card markers
    #     os.system(f"cp {cardmarkers_dir} {dep}")
    #     os.system(f"gzip -d ShortBRED_CARD_2017_markers.faa.gz")

    #download genomad database after downloading and activating the environment
    if filename == "genomad.yaml":
        env = os.path.join(outdir, "genomad")
        os.system(f"bash {genomad_switch} {conda_profile} {env}")