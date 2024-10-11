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


cat geno_pred/output/gwas_MAF250_piperacillin.tazobactam_bslmm_blup.param.txt



###
#concatenate results of GWAS
rm geno_pred/output/gwas_MAF250_ALL_bslmm_blup.param.txt
awk '{print $0, FILENAME}' geno_pred/output/gwas_MAF250_*_bslmm_blup.param.txt | awk '{gsub(/geno.*MAF250_|_bslmm.*/,"",$NF);print}' | sed '/^chr.*$/{x;/^$/!d;g;}'  > geno_pred/output/gwas_MAF250_ALL_bslmm_blup.param.txt

rm geno_pred/output/gwas_MAF250_ALL_bslmm_probit.param.txt
awk '{print $0, FILENAME}' geno_pred/output/gwas_MAF250_*_bslmm_probit.param.txt | awk '{gsub(/geno.*MAF250_|_bslmm.*/,"",$NF);print}' | sed '/^chr.*$/{x;/^$/!d;g;}'  > geno_pred/output/gwas_MAF250_ALL_bslmm_probit.param.txt

###
#concatenate results of GWAS presence absence
rm geno_pred/presence_absence/output/presence_absence_ALL_bslmm_blup.param.txt
awk '{print $0, FILENAME}' geno_pred/presence_absence/output/presence_absence_all_*_bslmm_blup.param.txt | awk '{gsub(/geno.*all_|_bslmm.*/,"",$NF);print}' | sed '/^chr.*$/{x;/^$/!d;g;}'  > geno_pred/presence_absence/output/presence_absence_ALL_bslmm_blup.param.txt
