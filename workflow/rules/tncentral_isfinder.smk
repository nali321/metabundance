OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]
TN_IS_DB = config["tn_is_db"]

rule tncentral_isfinder:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        tn_is=f"{OUTPUT}/tn_is/read_pair_{{sample}}/tn_is_{{sample}}.txt"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/blast
        blastn -query {input.assembly} \
        -db {TN_IS_DB} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out {OUTPUT}/tn_is/read_pair_{wildcards.sample}/tn_is_{wildcards.sample}.txt \
        '''