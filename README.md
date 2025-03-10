# Predicting antimicrobial resistance phenotypes across 7,000 E. coli genomes

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Context

At Arcadia, we're interested in mapping genotype-phenotype relationships at broader evolutionary scales than previously attempted. To do so, we're developing models that can capture genetic relationships—both linear and nonlinear—that may be inaccessible to conventional methods. A key part of this development is identifying and leveraging unique datasets that both capture large scales of diversity and rich phenotypic information.

In this repo, we characterize the genomic structure and perform genomic prediction using a diverse dataset of 7000 E. coli genomes we previously compiled and published. We first use exploratory population genomic analyses to verify that this dataset contains the type of high-quality genotypic information that can be leveraged for model development. We then follow this up with genomic prediction to uncover the genetic basis of three antimicrobial resistance (AMR) phenotypes. These genomic prediction analyses confirm and expand our understanding of the evolution of these AMR phenotypes and set the baseline for our future efforts, which will explore how non-linear models might be applied for similar genomic prediction goals.

## Data
Phenotypic and Genotypic data on the 7,000 E. coli genomes has previously been [published](https://research.arcadiascience.com/pub/dataset-ecoli-amr-genotype-phenotype/release/1#working-with-a-pangenome) and is available [here](https://zenodo.org/records/12692732)

Additional pre-computed data on presence-absence variation in the dataset is available [here](https://zenodo.org/records/14364732) and includes copies of files needed from the original 7,000 genome study needed to reprodouce all analyses.

## Installation and Setup
This repository uses Snakemake, R, and Python.
Dependency requirements are managed by conda.


The bash script `run_analyses.sh` can be used to initiate installation of miniforge3 (conda) and the main environment containing a snakemake installation, as well as running the snakemake pipeline to generate desired results. By default miniforge3 installation is commented out on the first lines of the script, if conda needs to be installed, uncomment these lines prior to running script/workflow. You may have to restart teminral after conda installation for conda to properly initialise for the first time.


## Basic workflow
The main snakemake file is found in `Ecoli_AMR_GenotypePhenotype/workflow/Snakemake`
This file will initialize downloading all necessary data, as well as trigerring all downstream analyses.
All dependencies are automiatically handled by conda using environments built from yaml files stored in Ecoli_AMR_GenotypePhenotype/workflow/envs


Analyses are split into 4 main subworkflows in Ecoli_AMR_GenotypePhenotype/workflow/rules
- `filtering.smk`:
This workflow cleans up the data, removing outlier samples, filtering down to informative sites
- `popgen_analyses.smk`:
This workflow generates some pop-gen visualizations such as a site-frequency spectrum and also constructs phylogenetic trees
- `genomic_prediction.smk`:
This workflow runs genomic prediction/GWAS using GEMMA
- `genomic_prediction_post_hoc.smk`:
  This workflow runs exploratory post-hoc genomic prediction analyses

## Directory Structure

### Supplemental pub tables
- `supplemental_tables`
Directory where supplemental pub tables can be accessed (both pdf and csv files available). Tables can also be reproduced by running snakemake pipeline (see above)

### Data
- `Ecoli_AMR_GenotypePhenotype/vcf_files`
Directory for storing genotypic information on SNPs/indels

- `Ecoli_AMR_GenotypePhenotype/presence_absence`
Directory for storing presence/absence data

- `Ecoli_AMR_GenotypePhenotype/pangenome`
Directory for storing data related to pangenome reference

- `Ecoli_AMR_GenotypePhenotype/phenotype_matrix`
Directory for storing data related to AMR phenotypes and assoaciated metadata


### Analyses

- `Ecoli_AMR_GenotypePhenotype/workflow/rules/scripts`
This is where most scripts are stored, generating output such as phylogenies, genomic prediction etc.

- `Ecoli_AMR_GenotypePhenotype/figs`
Directory where reproduced pub figures are saved

- `Ecoli_AMR_GenotypePhenotype/tables`
Directory where reproduced pub tables are saved

- `Ecoli_AMR_GenotypePhenotype/tree`
Directory for data related to phylogenetic analyses

- `Ecoli_AMR_GenotypePhenotype/geno_pred`
Directory for data and results from genomic prediction
