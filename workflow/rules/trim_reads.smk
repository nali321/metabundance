OUTPUT = config["output"]
ENVS = config["envs_path"]
ILLUMINACLIP = config["illuminaclip"]
CONDA_PATH = config["conda_path"]
ID = config["id"]
MUSCLE = config["muscle"]
USEARCH = config["usearch"]

rule trim_reads:
    input:
        forward_reads=config["forward"],
        reverse_reads=config["reverse"]
    output:
        forward_paired=f"{OUTPUT}/trimmed/forward_paired.fq.gz",
        forward_unpaired=f"{OUTPUT}/trimmed/forward_unpaired.fq.gz",
        reverse_paired=f"{OUTPUT}/trimmed/reverse_paired.fq.gz",
        reverse_unpaired=f"{OUTPUT}/trimmed/reverse_unpaired.fq.gz"
    shell:
        '''
        source {CONDA_PATH}
        conda activate {ENVS}/trimmomatic
        export JAVA_HOME={ENVS}/trimmomatic/bin/java
        trimmomatic PE -phred33 {input.forward_reads} {input.reverse_reads} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{ILLUMINACLIP}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
        conda deactivate
        '''