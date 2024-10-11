
grep -E 'gyr|parC|parE|bla|folA|Cluster_7017_|neo|ant1|group_4362'  pangenome/whole_pangenome.fasta | sed 's/>//g' | awk '{OFS='\t';print $1}' | less -S> geno_pred/gene_subset_bluphits.txt


####################################
#filter to specific genes
#bcftools view vcf_files/annotated_output.vcf.gz -r <(cat geno_pred/gene_subset.txt|  paste -s -d, -)   |less -S \
zgrep vcf_files/annotated_output.vcf.gz -f geno_pred/gene_subset_bluphits.txt -e '#'| \
bcftools view -s ^ERR1218638,ERR1218722,ERR4035496,ERR4037257,SRR10271534,SRR10272173,SRR10272272,SRR10272285,SRR10272401,SRR2449168,SRR3588897,SRR3932419,SRR3982215,SRR850803 | \
#bcftools filter -i 'AC > 7' | \
#bcftools +prune -m 0.25 -w 50000 | \
sed 's/\.\:\./0\:255\,0/g' | \

bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP[\t%GT]\n"  \
-o geno_pred/annotated_output_all_target_genes_tophits_geno.txt


####################################
#all common coding mutations
bcftools view vcf_files/annotated_output.vcf.gz -s ^ERR1218638,ERR1218722,ERR4035496,ERR4037257,SRR10271534,SRR10272173,SRR10272272,SRR10272285,SRR10272401,SRR2449168,SRR3588897,SRR3932419,SRR3982215,SRR850803 | \
bcftools filter -i 'AC > 300' | \
sed 's/\.\:\./0\:255\,0/g' > geno_pred/annotated_output_MAF300.vcf


#all common coding mutations
bcftools view vcf_files/annotated_output.vcf.gz -s ^ERR1218638,ERR1218722,ERR4035496,ERR4037257,SRR10271534,SRR10272173,SRR10272272,SRR10272285,SRR10272401,SRR2449168,SRR3588897,SRR3932419,SRR3982215,SRR850803 | \
bcftools filter -i 'AC > 379' | \
sed 's/\.\:\./0\:255\,0/g' | \
bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP[\t%GT]\n"  \
-o geno_pred/annotated_output_MAF380.txt




sed 's/\[[[:digit:]]\+\]\|\:GT//g' geno_pred/annotated_output_MAF380.txt | \
sed 's/\b0\b/-1/g' > geno_pred/annotated_output_MAF380_edited.txt

cat <(grep -e '#' geno_pred/annotated_output_biallelic_refedit_goodcontigs_ldpruned_remoutliers_edited.txt) <(grep -Ev '#' geno_pred/annotated_output_biallelic_refedit_goodcontigs_ldpruned_remoutliers_edited.txt | sed -n '0~10p') > geno_pred/annotated_output_biallelic_refedit_goodcontigs_ldpruned_remoutliers_edited_subsample.txt
