# metabundance

Metabundance is a pipeline designed to investigate metagenomic paired FASTQ reads at the assembly level and find abundance of antibiotic-resistance genes (ARGs) from within and between samples. It handles trimming the reads as well as creating and annotating the assemblies. The main feature of Metabundance is the ability for it to link each ARG it finds to various classes of mobile genetic elements (MGEs), as well as taxonomy.

Metabundance's output gives you two matrices of abundance data, as well as an observation matrix, to create [Phyloseq](https://joey711.github.io/phyloseq/) objects in R for further downstream analysis.

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
