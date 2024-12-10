OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
SNP_PATH = config["snp_path"]
SNP_DB = config["snp_db"]

rule kraken:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        taxa=f"{OUTPUT}/spraynpray/read_pair_{{sample}}/read_pair_{{sample}}.csv"
    shell:
        '''
        set +eu
        source {CONDA_PATH}
        conda activate {SNP_PATH}
        spray-and-pray.py -g {input.assembly} \
        -out {OUTPUT}/spraynpray/read_pair_{wildcards.sample} \
        -ref {SNP_DB}
        conda deactivate
        set -eu
        '''