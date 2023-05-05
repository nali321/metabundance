# metabundance

Metabundance is a pipeline utilizing Snakemake designed to investigate metagenomic paired FASTQ reads at the assembly level and find abundance of antibiotic-resistance genes (ARGs) from within and between samples. It handles trimming the reads as well as creating and annotating the assemblies. The main feature of Metabundance is the ability for it to link each ARG it finds to various classes of mobile genetic elements (MGEs), as well as taxonomy.

Metabundance's output gives you two matrices of abundance data, as well as an observation matrix, to create [Phyloseq](https://joey711.github.io/phyloseq/) objects in R for further downstream analysis and visualizations.

## Prerequisites:
- Must be ran on Linux
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed and conda set up
- [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) installed
- Read pairs need to be labelled properly:
  - numbered 1 - sample total following this format: [sample number]_[R1/R2]_001.fastq.gz
- If you want taxonomy data, you need to download [SprayNPray](https://github.com/Arkadiy-Garber/SprayNPray) as a conda environment

## Installation:
```
git clone https://github.com/nali321/metabundance
```

## Setup:
Run conda_installer.py to create environments for pipeline.

```
python /path/to/conda_installer.py -o /path/to/envs
conda activate /path/to/envs/metabundance
```

This pipeline was originally designed on an HPC-server without internet access on clusters, therefore pre-installing the conda environments/dependencies before running the Snakemake pipeline was implemented. This achieves the same result as if the conda module was used on Snakemake, where the conda environments are specifically downloaded to a separate folder instead of your main envs folder.

## Usage Guide:
### Quick Start
```
python /path/to/metabundance.py -r /path/to/paired_reads -e /path/to/envs -i /path/to/illuminaclip -m /path/to/metadata.csv -o /path/to/outdir
```
