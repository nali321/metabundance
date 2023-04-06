OUTPUT = config["output"]
ENVS = config["envs_path"]
CONDA_PATH = config["conda_path"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]
DIR = config["reads"]

rule shortbred_quantify:
    input:
        markers=f"{OUTPUT}/shortbred_identify/read_pair_{{sample}}/markers.faa",
        forward_reads=f"{DIR}/{{sample}}_R1_001.fastq.gz",
        reverse_reads=f"{DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        folder=f"{OUTPUT}/shortbred_quantify/read_pair_{{sample}}/tmp_quantify",
        abundance=f"{OUTPUT}/shortbred_quantify/read_pair_{{sample}}/results.tsv"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/shortbred
        shortbred_quantify.py --markers {input.markers} --wgs {input.forward_reads} \
        {input.reverse_reads} --results {output.abundance} \
        --tmp {output.folder} \
        --usearch {USEARCH}
        conda deactivate
        '''