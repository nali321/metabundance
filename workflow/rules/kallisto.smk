OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]

rule kallisto:
    input:
        idx=f"{OUTPUT}/kallisto/args.idx",
        forward_reads=f"{DIR}/read_pair_{{sample}}_R1_001.fastq.gz",
        reverse_reads=f"{DIR}/read_pair_{{sample}}_R2_001.fastq.gz"
    output:
        abundance=f"{OUTPUT}/kallisto/read_pair_{{sample}}/abundance.tsv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/kallisto
        kallisto quant -i {input.idx} -o {OUTPUT}/kallisto/read_pair_{wildcards.sample} {input.forward_reads} {input.reverse_reads}
        conda deactivate
        '''