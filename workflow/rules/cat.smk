OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
CAT_PROFILE = config["cat_profile"]
CAT_PATH = config["cat_path"]
CAT_PACK = config["cat_pack"]
CAT_DB = config["cat_db"]

rule cat:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        taxa=f"{OUTPUT}/cat/read_pair_{{sample}}/CAT.contigs2classification.txt"
    shell:
        '''
        source {CAT_PROFILE}
        conda activate {CAT_PATH}
        {CAT_PACK} contigs \
        -c {input.assembly} \
        -d {CAT_DB}/db \
        -t {CAT_DB}/tax  \
        -n 16 -o {OUTPUT}/cat/read_pair_{wildcards.sample} \
        --out_prefix {OUTPUT}/cat/read_pair_{wildcards.sample}/CAT
        conda deactivate
        '''