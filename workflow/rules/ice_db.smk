OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]
ICE_DB = config["ice_db"]

rule ice_db:
    input:
        db=f"{ICE_DB}"
    output:
        ice_db=f"{OUTPUT}/database/ice/ICE_seq_all.fas"
    shell:
        '''
        mkdir -p {OUTPUT}/ice/db
        cp {ICE_DB} {OUTPUT}/ice/db
        source {CONDA_PATH}
        conda activate {ENVS}/blast_env
        makeblastdb -in {OUTPUT}/ice/db/ICE_seq_all.fas -dbtype nucl
        conda deactivate
        '''
