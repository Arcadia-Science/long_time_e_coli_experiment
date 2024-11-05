#install miniforge3 and restart terminal
#curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

#main wd of projectm make output directories for pipeline to run
cd long_time_e_coli_experiment/Ecoli_AMR_GenotypePhenotype
mkdir figs
mkdir tables

#wd of snakemake pipeline
cd l/workflow

#setup and activate conda environment for snakemake to run
conda env create --file=envs/snakemake.yaml
conda activate snakemake

#dry run test
snakemake -n
#build dag of pipeline
snakemake --forceall --dag | dot -Tpdf > dag.pdf
