# metabundance

Metabundance is a pipeline designed to investigate metagenomic paired FASTQ reads at the assembly level and find abundance of antibiotic-resistance genes (ARGs) from within and between samples. It handles trimming the reads as well as creating and annotating the assemblies.

Metabundance output is supposed to be used primarily with Phyloseq in R. It gives you two matrices (samples as columns, genes as rows) of abundance data, as well as an observation matrix (gene information as columns, genes as rows), that are to be used as the inputs to create Phyloseq objects. 

Metabundance also functions as a way to link together ARG contigs with mobile genetic elements (MGEs), as well as finding the taxonomy of each contig that an ARG is on.

Metabundance uses many depenencies in its pipeline, therefore a script (conda_installer.py) is used to download every single dependency and save it within a separate file containing conda environments (much like how Snakemake can create temporary conda environments).

## Prerequisites:
- Must be ran on Linux
- Must have Miniconda installed
- Must have Mamba installed
- Read pairs need to be labelled properly:
  - numbered 1 - sample total following this format: [sample number]_[R1/R2]_001.fq.gz
- If you want taxonomy data, you need to download [SprayNPray](https://github.com/Arkadiy-Garber/SprayNPray) as a conda environment

## Usage Guide:
### Setup
1. Run conda_installer.py to create environments for pipeline

```
python /path/to/conda_installer.py stuff
```

2. Run reads2args.py to create metagenomic assemblies and annotate them

3. Run args2abundance.py to find ARG abundance and create matrices needed for Phyloseq
