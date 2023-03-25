OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

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