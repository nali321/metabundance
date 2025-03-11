OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
CAT_PATH = config["cat_path"]
CAT_DB = config["cat_db"]
CAT_CONDA = config["cat_conda"]

rule cat:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        taxa=f"{OUTPUT}/cat/read_pair_{{sample}}/CAT_{{sample}}.contigs2classification.txt"
    shell:
        '''
        source {CAT_CONDA}
        conda activate {CAT_PATH}
        {CAT_PATH} contigs \
        -c {input.assembly} \
        -d {CAT_DB}/db \
        -t {CAT_DB}/tax  \
        -n 16 -o {OUTPUT}/cat/read_pair_{wildcards.sample} \
        --out_prefix {OUTPUT}/cat/read_pair_{wildcards.sample}/CAT_{wildcards.sample}
        conda deactivate
        '''