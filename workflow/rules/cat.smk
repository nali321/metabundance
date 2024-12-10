OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
CAT_PATH = config["cat_path"]
CAT_DB = config["cat_db"]

rule cat:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        taxa=f"{OUTPUT}/cat/read_pair_{{sample}}/read_pair_{{sample}}.csv"
    shell:
        '''
        source {CONDA_PATH}
        {CAT_PATH} contigs \
        -c {input.assembly} \
        -d {cat_db}/db \
        -t {cat_db}/tax  \
        -n 16 -o {OUTPUT}/cat/read_pair_{wildcards.sample} \
        --out_prefix {OUTPUT}/cat/read_pair_{wildcards.sample}/CAT_{wildcards.sample}
        conda deactivate
        '''