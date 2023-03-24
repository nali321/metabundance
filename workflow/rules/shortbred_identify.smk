OUTPUT = config["output"]
ENVS = config["envs_path"]
ILLUMINACLIP = config["illuminaclip"]
CONDA_PATH = config["conda_path"]
ID = config["id"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]

rule shortbred_identify:
    input:
        proteins=config["faa_file"],
        ref_markers=config["CARD_markers"],
    output:
        markers=f"{OUTPUT}/shortbred/shortbred_identify/read_pair_{ID}/markers.faa"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/shortbred
        shortbred_identify.py --goi {input.proteins} \
        --ref {input.ref_markers} --tmp {OUTPUT}/shortbred/shortbred_identify/read_pair_{ID} \
        --markers {output.markers} --usearch {USEARCH} \
        --muscle {MUSCLE} \
        --cdhit {ENVS}/cdhit/bin/cd-hit \
        --blastp {ENVS}/blast_env/bin/blastp \
        --makeblastdb {ENVS}/blast_env/bin/makeblastdb
        conda deactivate
        '''