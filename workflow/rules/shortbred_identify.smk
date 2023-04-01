OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]
FASTA = config["fasta"]

rule shortbred_identify:
    input:
        proteins=f"{FASTA}/protein_files/{{sample}}.faa",
        ref_markers=config["CARD_markers"],
    output:
        markers=f"{OUTPUT}/shortbred_identify/read_pair_{{sample}}/markers.faa"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/shortbred
        shortbred_identify.py --goi {input.proteins} \
        --ref {input.ref_markers} --tmp {OUTPUT}/shortbred_identify/read_pair_{wildcards.sample} \
        --markers {output.markers} --usearch {USEARCH} \
        --muscle {MUSCLE} \
        --cdhit {ENVS}/cdhit/bin/cd-hit \
        --blastp {ENVS}/blast_env/bin/blastp \
        --makeblastdb {ENVS}/blast_env/bin/makeblastdb
        conda deactivate
        '''