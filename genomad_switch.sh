#!/bin/bash

source $1

conda activate $2

cd $3

genomad download-database .

conda deactivate