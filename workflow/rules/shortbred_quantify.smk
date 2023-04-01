OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]

rule shortbred_quantify:
    input:
        markers=f"{OUTPUT}/shortbred/shortbred_identify/read_pair_{ID}/markers.faa",
        forward_reads=config["forward"],
        reverse_reads=config["reverse"]
    output:
        abundance=f"{OUTPUT}/shortbred/shortbred_quantify/read_pair_{ID}/results.tsv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/shortbred
        shortbred_quantify.py --markers {input.markers} --wgs {input.forward_reads} \
        {input.reverse_reads} --results {output.abundance} \
        --tmp {OUTPUT}/shortbred/shortbred_quantify/read_pair_{ID}/tmp_quantify \
        --usearch {USEARCH}
        conda deactivate
        '''