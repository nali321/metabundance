OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]
ICE_DB = config["ice_db"]

rule ice_db:
    input:
        db=ICE_DB
    output:
        ice_db=f"{OUTPUT}/db/ice/ICE_seq_all.fas"
    shell:
        '''
        mkdir -p {OUTPUT}/db/ice
        cp {input.db} {OUTPUT}/db/ice
        source {CONDA_PATH}
        conda activate {ENVS}/blast_env
        makeblastdb -in {OUTPUT}/db/ice/ICE_seq_all.fas -dbtype nucl
        conda deactivate
        '''
