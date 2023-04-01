OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule genomad:
    input:
        assembly=f"{OUTPUT}/metaspades/{{sample}}/contigs.fasta"
    output:
        plasmid=f"{OUTPUT}/genomad/{{sample}}/contigs_summary/contigs_plasmid_summary.tsv",
        virus=f"{OUTPUT}/genomad/{{sample}}/contigs_summary/contigs_virus_summary.tsv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/genomad
        genomad end-to-end {input.assembly} {OUTPUT}/genomad/{wildcards.sample} {ENVS}/genomad/genomad_db
        '''