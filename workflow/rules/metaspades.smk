OUTPUT = config["output"]
ENVS = config["envs_path"]
ILLUMINACLIP = config["illuminaclip"]
CONDA_PATH = config["conda_path"]
ID = config["id"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]

rule metaspades:
    input:
        one=f"{OUTPUT}/trimmed/forward_paired.fq.gz",
        two=f"{OUTPUT}/trimmed/reverse_paired.fq.gz",
        s=f"{OUTPUT}/trimmed/forward_unpaired.fq.gz"
    output:
        assembly=f"{OUTPUT}/metaspades/contigs.fasta"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/spades
        metaspades.py -o {OUTPUT}/metaspades -1 {input.one} -2 {input.two} -s {input.s} -m 128 -k 33,55,77,99 -t 16
        conda deactivate
        '''