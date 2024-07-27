#!/bin/bash

cd Ecoli_AMR_GenotypePhenotype/



#unannotated files
bcftools filter -e 'TYPE = "indel"' vcf_files/filtered_output.vcf.gz |  \
bcftools view -M2 -m2 -v snps -Oz  > vcf_files/filtered_output_biallelic.vcf.gz

bcftools query vcf_files/filtered_output_biallelic.vcf.gz  -f "%CHROM\t%POS\t%ALT\t[\t%GT]\n" \
-o vcf_files/filtered_output_biallelic.txt





#########################################################
#annotated files


bcftools filter -e 'TYPE = "indel"' vcf_files/annotated_output.vcf.gz |  \
bcftools view -M2 -m2 -v snps -Oz  > vcf_files/annotated_output_biallelic.vcf.gz

#bcftools query vcf_files/annotated_output_biallelic.vcf.gz  -f "%CHROM\t%POS\t%ALT\t[\t%GT]\n" \
#-o vcf_files/annotated_output_biallelic.txt

#############################
#synonymous bi allelic sites only
zgrep -E 'synonymous_variant|#' vcf_files/annotated_output_biallelic.vcf.gz > vcf_files/annotated_output_biallelic_synonymous.vcf

#############
bcftools filter -i 'AC > 7' vcf_files/annotated_output_biallelic_synonymous.vcf | sed 's/\.\:\./0\:255\,0/g'  > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit.vcf

cat <(grep -e '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit.vcf) <(grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit.vcf | sed -n '0~10p') > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_subsample.vcf


#grep -e '#' vcf_files/annotated_output_biallelic_synonymous.vcf > vcf_files/annotated_output_biallelic_synonymous_header.vcf
#grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous.vcf | sed -n '0~50p' > vcf_files/annotated_output_biallelic_synonymous_subsample_nohead.vcf
#cat vcf_files/annotated_output_biallelic_synonymous_header.vcf vcf_files/annotated_output_biallelic_synonymous_subsample_nohead.vcf >  vcf_files/annotated_output_biallelic_synonymous_subsample.vcf

rm vcf_files/filtered_output_biallelic_header.vcf
rm vcf_files/subsample_nohead.vcf


##################################################
#SFS data prep
zgrep -E 'synonymous|missense|HIGH|#' vcf_files/annotated_output_biallelic.vcf.gz | \
bcftools query   -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t%ANN\n" -o vcf_files/annotated_output_biallelic_SFS_counts.txt

#sum value of genotype fields to get allele count
#awk '{for(i=4;i<=NF;i++) t+=$i; print $1,$2,$3,t,$NF; t=0}' vcf_files/annotated_output_biallelic_SFS.txt > vcf_files/annotated_output_biallelic_SFS_counts.txt
