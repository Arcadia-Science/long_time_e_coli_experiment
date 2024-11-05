#!/bin/bash

#get data

aws s3 cp s3://ecoli-pangenomics-data/Ecoli_AMR_GenotypePhenotype/ Ecoli_AMR_GenotypePhenotype --recursive

#importing updates metadata
#aws s3 cp s3://ecoli-pangenomics-data/Ecoli_AMR_GenotypePhenotype/strains_metadata.csv Ecoli_AMR_GenotypePhenotype/strains_metadata.csv
#aws s3 cp s3://ecoli-pangenomics-data/Ecoli_AMR_GenotypePhenotype/phenotype_matrix Ecoli_AMR_GenotypePhenotype/phenotype_matrix --recursive
