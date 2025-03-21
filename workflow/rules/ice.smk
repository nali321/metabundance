OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]
ICE_DB = config["ice_db"]

rule ice:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        ice=f"{OUTPUT}/ice/read_pair_{{sample}}/ice.txt"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/blast_env
        blastn -query {input.assembly} \
        -db {ICE_DB} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out {OUTPUT}/ice/read_pair_{wildcards.sample}/ice.txt
        '''