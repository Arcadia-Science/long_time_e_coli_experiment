#!/bin/bash

#create binary ped file for input to gemma using plink
plink --vcf geno_pred/annotated_output_MAF300.vcf \
--make-bed --out geno_pred/plink/annotated_output_MAF300 --allow-extra-chr


# gemma calculate relatedness matrix
##1 is centered relmatrix better for pop structure in 'lower organisms'
##2 is standardized relmatrix preferred if rarer SNPs expected to have larger effects
gemma -bfile geno_pred/plink/annotated_output_MAF300  -gk 1 -o annotated_output_MAF300





#naive linear model with no relatedness matrix
#gemma -bfile geno_pred/plink/annotated_output_MAF500 -lm 2 -o gwas_MAF500

#linear mixed model with relatedness matrix LRT
gemma -bfile geno_pred/plink/annotated_output_MAF500 -lmm 2 -k output/annotated_output_MAF500.cXX.txt -o gwas_MAF500

#linear mixed model with relatedness matrix LRT
gemma -bfile geno_pred/plink/annotated_output_MAF500 -lmm 1 -k output/annotated_output_MAF500.cXX.txt -o gwas_MAF500_wald





#bayesian sparse linear mixed models
#3 is probit model
gemma -bfile geno_pred/plink/annotated_output_MAF300 -bslmm 3 -k output/annotated_output_MAF300.cXX.txt -o gwas_MAF300_ampicillin_bslmm_probit

#2 is ridge/BLUP fastest
gemma -bfile geno_pred/plink/annotated_output_MAF300 -bslmm 2 -k output/annotated_output_MAF300.cXX.txt -o gwas_MAF300_ampicillin_bslmm_blup
