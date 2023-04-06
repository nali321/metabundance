OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
FASTA = config["fasta"]

rule kallisto_index:
    input:
        args=FASTA
    output:
        idx=f"{OUTPUT}/kallisto/args.idx"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/kallisto
        kallisto index -i {output.idx} {input.args}
        conda deactivate
        '''