OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule integron_finder:
    input:
        assembly=f"{OUTPUT}/metaspades/read_pair_{{sample}}/contigs.fasta"
    output:
        integrons=f"{OUTPUT}/integron_finder/read_pair_{{sample}}/Results_Integron_Finder_contigs/contigs.integrons"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/integron_finder
        integron_finder {input.assembly} --outdir {OUTPUT}/integron_finder/read_pair_{wildcards.sample}
        '''