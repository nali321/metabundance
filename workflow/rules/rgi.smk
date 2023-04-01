OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
ID = config["id"]

rule rgi:
    input:
        assembly=f"{OUTPUT}/metaspades/{{sample}}/contigs.fasta"
    output:
        args=f"{OUTPUT}/rgi/{{sample}}/final.txt"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/rgi
        cd {OUTPUT}/rgi/{{sample}}
        rgi main --input_sequence {input.assembly} --output_file final -a DIAMOND --input_type contig --include_loose
        conda deactivate
        '''