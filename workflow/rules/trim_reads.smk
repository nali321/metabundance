OUTPUT = config["output"]
ENVS = config["envs_path"]
ILLUMINACLIP = config["illuminaclip"]
CONDA_PATH = config["conda_path"]
DIR = config["reads"]

#export JAVA_HOME={ENVS}/trimmomatic/bin/java

rule trim_reads:
    input:
        forward_reads=f"{DIR}/{{sample}}_R1_001.fastq.gz",
        reverse_reads=f"{DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        forward_paired=f"{OUTPUT}/trimmed/read_pair_{{sample}}/forward_paired.fq.gz",
        forward_unpaired=f"{OUTPUT}/trimmed/read_pair_{{sample}}/forward_unpaired.fq.gz",
        reverse_paired=f"{OUTPUT}/trimmed/read_pair_{{sample}}/reverse_paired.fq.gz",
        reverse_unpaired=f"{OUTPUT}/trimmed/read_pair_{{sample}}/reverse_unpaired.fq.gz"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/trimmomatic
        export JAVA_HOME={ENVS}/trimmomatic
        trimmomatic PE -phred33 {input.forward_reads} {input.reverse_reads} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{ILLUMINACLIP}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
        conda deactivate
        '''