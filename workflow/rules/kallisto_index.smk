OUTPUT = config["output"]
ENVS = config["envs_path"]
ILLUMINACLIP = config["illuminaclip"]
CONDA_PATH = config["conda_path"]
ID = config["id"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]

rule kallisto_index:
    input:
        args=config["all_args"]
    output:
        idx=f"{OUTPUT}/kallisto/args.idx"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/kallisto
        kallisto index -i {output.idx} {input.args}
        conda deactivate
        '''