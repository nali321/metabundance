#!/bin/bash

source $1

conda activate $2

genomad download-database $2

conda deactivate