#!/bin/bash

#get data
#TODO change to public zenodo
aws s3 cp s3://ecoli-pangenomics-data/Ecoli_AMR_GenotypePhenotype/ Ecoli_AMR_GenotypePhenotype --recursive

