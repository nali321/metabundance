OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule genomad:
    input:
        assembly=f"{OUTPUT}/metaspades/{{sample}}/contigs.fasta"
    output:
        plasmid=f"{OUTPUT}/rgi/final.txt"
        virus=f"{OUTPUT}/genomad/{{sample}}/"
    shell:
        '''
        source {CONDA_PATH}
        conda activate genomad
        cp {input.assembly} ~/wastewater_samples_1/concat/genomad/read_pair_"$i"
        cd ~/wastewater_samples_1/concat/genomad/read_pair_"$i"
        genomad end-to-end contigs.fasta {{sample}} /mmfs1/home/4565alin/build/genomad/genomad_db
        '''