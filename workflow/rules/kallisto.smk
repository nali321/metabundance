OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
ID = config["id"]

rule kallisto:
    input:
        idx=f"{OUTPUT}/kallisto/args.idx",
        forward_reads=config["forward"],
        reverse_reads=config["reverse"]
    output:
        abundance=f"{OUTPUT}/kallisto/read_pair_{ID}/abundance.tsv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/kallisto
        kallisto quant -i {input.idx} -o {OUTPUT}/kallisto/read_pair_{ID} {input.forward_reads} {input.reverse_reads}
        conda deactivate
        '''