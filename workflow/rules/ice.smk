OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]
ICE_DB = config["ice_db"]

rule ice:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta",
        db=f"{OUTPUT}/ice/db/ICE_seq_all.fas"
    output:
        ice=f"{OUTPUT}/ice/read_pair_{{sample}}/ice.txt"
    shell:
        '''
        mkdir -p {OUTPUT}/ice/read_pair_{wildcards.sample}
        source {CONDA_PATH}
        conda activate {ENVS}/blast_env
        blastn -query {input.assembly} \
        -db {input.db} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out {OUTPUT}/ice/read_pair_{wildcards.sample}/ice.txt
        conda deactivate
        '''