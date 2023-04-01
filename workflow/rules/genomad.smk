OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule genomad:
    input:
        assembly=f"{OUTPUT}/metaspades/contigs.fasta"
    output:
        args=f"{OUTPUT}/rgi/final.txt"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/rgi
        cd {OUTPUT}/rgi
        rgi main --input_sequence {input.assembly} --output_file final -a DIAMOND --input_type contig --include_loose
        conda deactivate
        '''