#!/bin/bash

#calculate r2 ld statistic
plink --vcf vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_subsample.vcf \
--make-bed --out plink/annotated_output_biallelic_synonymous_MAF7_refedit_subsample \
--r2 --allow-extra-chr --ld-window-r2 0 --ld-window-kb 100000


#sample ped file
plink --vcf tree/annotated_output_bla_refedit.vcf \
--recode --out TEST_PLINK --allow-extra-chr 