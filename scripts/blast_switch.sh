#!/bin/bash

source $1

mamba create -p $2

conda activate $2

conda install -y -c bioconda blast

conda deactivate