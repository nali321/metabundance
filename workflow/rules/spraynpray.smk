OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]

rule spraynpray:
    input:
        assembly=f"{OUTPUT}/metaspades/{{sample}}/contigs.fasta"
    output:
        taxa=f"{OUTPUT}/spraynpray/{{sample}}/{{sample}}.csv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/sprayandpray
        spray-and-pray.py -g {assembly.input} \
        -out {OUTPUT}/spraynpray/{wildcards.sample} \
        -ref /mmfs1/home/4565alin/build/spraynpraydb/nr.faa
        ''''