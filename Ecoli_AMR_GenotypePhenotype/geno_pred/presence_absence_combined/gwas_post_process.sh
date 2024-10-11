#concatenate results of GWAS presence absence full BLUP
rm geno_pred/presence_absence_all/output/presence_absence_ALL_bslmm_blup.param.txt
awk '{print $0, FILENAME}' geno_pred/presence_absence_combined/output/annotated_output_MAF250_presence_absence_all_*_bslmm_blup.param.txt| awk '{gsub(/geno.*all_|_bslmm.*/,"",$NF);print}' | sed '/^chr.*$/{x;/^$/!d;g;}'  > geno_pred/presence_absence_combined/output/presence_absence_ALL_bslmm_blup.param.txt

#same but for probit model on focal 4 phenotypes
rm geno_pred/presence_absence_all/output/presence_absence_ALL_bslmm_probit.param.txt
awk '{print $0, FILENAME}' geno_pred/presence_absence_combined/output/annotated_output_MAF250_presence_absence_all_*_bslmm_probit.param.txt| awk '{gsub(/geno.*all_|_bslmm.*/,"",$NF);print}' | sed '/^chr.*$/{x;/^$/!d;g;}'  > geno_pred/presence_absence_combined/output/presence_absence_ALL_bslmm_blup.param.txt
