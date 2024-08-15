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

#only good pangenome contigs
cat <(grep -e '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit.vcf) \
 <(grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit.vcf  | \
 grep -f pangenome/ecor_shared_contigs.txt ) > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs.vcf

cat <(grep -e '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs.vcf) <(grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs.vcf | sed -n '0~10p') > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.vcf


bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP[\t%GT]\n" \
vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.vcf \
-o vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample_geno.txt



#######################
#remove outlier samples
cat <(zgrep -e '#' vcf_files/annotated_output_biallelic.vcf.gz) \
 <(zgrep -Ev '#' vcf_files/annotated_output_biallelic.vcf.gz | zgrep -f pangenome/ecor_shared_contigs.txt ) | \
zgrep -E 'synonymous|missense|HIGH|#' | \
bcftools view -Ov -s ^ERR1218638,ERR1218722,ERR4035496,ERR4037257,SRR10271534,SRR10272173,SRR10272272,SRR10272285,SRR10272401,SRR2449168,SRR3588897,SRR3932419,SRR3982215,SRR850803 \
-o vcf_files/annotated_output_biallelic_goodcontigs_remoutliers.vcf


##################################################
#SFS data prep
zgrep -E 'synonymous|missense|HIGH|#' vcf_files/annotated_output_biallelic.vcf.gz | \
bcftools query   -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t%ANN\n" -o vcf_files/annotated_output_biallelic_SFS_counts.txt




#SFS data prep good contigs and remove quad outlier samples
cat <(zgrep -e '#' vcf_files/annotated_output_biallelic.vcf.gz) \
 <(zgrep -Ev '#' vcf_files/annotated_output_biallelic.vcf.gz | zgrep -f pangenome/ecor_shared_contigs.txt ) | \
zgrep -E 'synonymous|missense|HIGH|#' | \
bcftools view -s ^ERR1218638,ERR1218722,ERR4035496,ERR4037257,SRR10271534,SRR10272173,SRR10272272,SRR10272285,SRR10272401,SRR2449168,SRR3588897,SRR3932419,SRR3982215,SRR850803 | \
bcftools query  -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t%ANN\n" -o vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts_remoutliers.txt


#sum value of genotype fields to get allele count
#awk '{for(i=4;i<=NF;i++) t+=$i; print $1,$2,$3,t,$NF; t=0}' vcf_files/annotated_output_biallelic_SFS.txt > vcf_files/annotated_output_biallelic_SFS_counts.txt


#looking into quadrupletons
cat <(zgrep -e '#' vcf_files/annotated_output_biallelic.vcf.gz) \
 <(zgrep -Ev '#' vcf_files/annotated_output_biallelic.vcf.gz | zgrep -f pangenome/ecor_shared_contigs.txt ) | \
zgrep -E 'synonymous|#' | \
bcftools filter -i 'INFO/AC = 4'| sed 's/\.\:\./0\:255\,0/g'  > vcf_files/annotated_output_biallelic_synonymous_AC4_refedit_goodcontigs.vcf


cat <(zgrep -e '#' vcf_files/annotated_output_biallelic.vcf.gz) \
 <(zgrep -Ev '#' vcf_files/annotated_output_biallelic.vcf.gz | zgrep -f pangenome/ecor_shared_contigs.txt ) | \
zgrep -E 'synonymous|#' | \
bcftools filter -i 'INFO/AC = 4'| \
bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t[\t%GT]\t%ANN\n" -o vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt





###############################
#LD pruning

bcftools +prune -m 0.25 -w 50000 vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs.vcf \
-Ov > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_ldpruned.vcf


cat <(grep -e '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_ldpruned.vcf) <(grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_ldpruned.vcf | sed -n '0~10p') > vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_ldpruned_subsample.vcf
