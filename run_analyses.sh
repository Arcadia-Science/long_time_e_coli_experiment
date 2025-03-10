#!/bin/bash

set -e

#install miniforge3 and restart terminal
#curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

#main wd of projectm make output directories for pipeline to run
cd Ecoli_AMR_GenotypePhenotype
mkdir figs
mkdir tables

#wd of snakemake pipeline
cd workflow

#setup and activate conda environment for snakemake to run
conda env create --file=envs/snakemake.yaml
conda activate snakemake

#dry run test
snakemake -n
#build dag of pipeline
snakemake --forceall --dag | dot -Tpdf > dag.pdf
#run snakemake pipeline
snakemake --use-conda
