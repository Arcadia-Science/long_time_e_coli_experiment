plink --bfile geno_pred/presence_absence_combined/plink/annotated_output_MAF250_presence_absence_all  --make-bed --allow-extra-chr \
    --extract range  gwas_presence_absence_complete_bslmm_probit_top10_ranges.txt \
    --out geno_pred/plink/presence_absence_complete_probit_tophits

#sed -i 's/-9/2/g' geno_pred/plink/presence_absence_complete_probit_tophits.fam


awk '$2=$4' geno_pred/plink/TESTBED.bim > geno_pred/plink/TESTBED_tmp.bim
mv geno_pred/plink/TESTBED_tmp.bim geno_pred/plink/TESTBED.bim


cp geno_pred/plink/annotated_output_MAF250_ciprofloxacin.fam geno_pred/plink/TESTBED.fam
#rm geno_pred/plink/TESTBED_tmp.bim

awk  '{gsub("1","2",$6)}1' geno_pred/plink/annotated_output_MAF250_ciprofloxacin.fam |awk  '{gsub("0","1",$6)}1' >geno_pred/plink/TESTBED.fam



plink --bfile geno_pred/plink/TESTBED --epistasis --allow-extra-chr --allow-no-sex --out geno_pred/plink/TESTBED --epi1 1




##############


bcftools filter geno_pred/annotated_output_MAF250.vcf -T geno_pred/gwas_MAF250_ALL_bslmm_probit_top_hits_sites.txt | \
bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP[\t%GT]\n"  \
-o geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt



bcftools filter geno_pred/annotated_output_MAF250.vcf -T geno_pred/gwas_MAF250_ALL_bslmm_probit_top_hits_sites.txt > geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits.vcf


plink --bfile geno_pred/plink/presence_absence_complete_probit_tophits --r2 --allow-extra-chr --ld-window-r2 0 --inter-chr --out geno_pred/plink/presence_absence_complete_probit_tophits

#null distribution of SNPs

#sample 3k snpsindels first
awk '{print $2}' geno_pred/presence_absence_combined/plink/annotated_output_MAF250_presence_absence_all.bim  | grep -Ev 'pres' | \
shuf -n 3000 > geno_pred/presence_absence_combined_ld_marker_null_set.txt
#sample 1k presabs markers next
awk '{print $2}' geno_pred/presence_absence_combined/plink/annotated_output_MAF250_presence_absence_all.bim  | grep 'pres' | \
shuf -n 1000 >> geno_pred/presence_absence_combined_ld_marker_null_set.txt


#extract null ld SNPs
plink --bfile geno_pred/presence_absence_combined/plink/annotated_output_MAF250_presence_absence_all \
    --make-bed --extract geno_pred/presence_absence_combined_ld_marker_null_set.txt \
    --out geno_pred/presence_absence_combined_ld_marker_null_set \
    --allow-extra-chr

plink --bfile geno_pred/presence_absence_combined_ld_marker_null_set --r2 --allow-extra-chr --ld-window-r2 0 --inter-chr --out geno_pred/presence_absence_combined_ld_marker_null_set
