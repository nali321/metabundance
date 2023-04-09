OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule metaspades:
    input:
        one=f"{OUTPUT}/trimmed/read_pair_{{sample}}/forward_paired.fq.gz",
        two=f"{OUTPUT}/trimmed/read_pair_{{sample}}/reverse_paired.fq.gz",
        s=f"{OUTPUT}/trimmed/read_pair_{{sample}}/forward_unpaired.fq.gz"
    output:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/spades
        metaspades.py -o {OUTPUT}/metaspades/read_pair_{wildcards.sample} -1 {input.one} -2 {input.two} -s {input.s} -m 128 -k 33,55,77,99 -t 16
        conda deactivate
        '''