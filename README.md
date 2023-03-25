# metabundance

Metabundance is a pipeline designed to trim, assemble, and annotate metagenomic assemblies from paired fastq files, then find abundance of antibiotic-resistance genes (ARGs) from within and between samples.

Metabundance output is supposed to be used primarily with Phyloseq in R. It gives you two matrices (samples as columns, genes as rows) of abundance data, as well as an observation matrix (gene information as columns, genes as rows), that are to be used as the inputs to create Phyloseq objects. 

Metabundance also functions as a way to link together ARG contigs with mobile genetic elements (MGEs), as well as finding the taxonomy of each contig that an ARG is on.

Metabundance uses many depenencies in its pipeline, therefore a script (conda_installer.py) is used to download every single dependency and save it within a separate file containing conda environments (much like how Snakemake can create temporary conda environments).

PREQUISITIES:
- Must be ran on Linux
- Must have Miniconda installed
- Must have Mamba installed

STEPS:

1. Run conda_installer.py to create environments for pipeline

2. Run reads2args.py to create metagenomic assemblies and annotate them

3. Run args2abundance.py to find ARG abundance and create matrices needed for Phyloseq
