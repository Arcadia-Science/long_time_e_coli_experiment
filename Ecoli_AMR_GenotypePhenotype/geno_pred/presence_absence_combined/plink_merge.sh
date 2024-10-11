
#first prep files for merging by giving unique SNP ids

#presence absence files
#for presence absence data, first give bogus position, then unique ID based on contig name and presabs in bim file
awk '{gsub(/.*/,"696969",$4);print}' geno_pred/presence_absence/presence_absence_all.bim  | \
awk '$2=$1 "presabs"' > geno_pred/presence_absence_combined/plink/presence_absence_all.bim

#move over bed file
cp geno_pred/presence_absence/presence_absence_all.bed geno_pred/presence_absence_combined/plink/presence_absence_all.bed
cp geno_pred/presence_absence/presence_absence_all.fam geno_pred/presence_absence_combined/plink/presence_absence_all.fam

#similary edit map file
awk '{gsub(/.*/,"69696969",$3);print}' geno_pred/presence_absence/presence_absence_all.map | \
awk '$2=$1 "presabs"' > geno_pred/presence_absence_combined/plink/presence_absence_all.map


#SNP and indel files
#edit bim file to SNP id is contracted contig name, plus line number and position to make it unique
awk '$2=$1' geno_pred/plink/annotated_output_MAF250.bim |awk '{gsub(/_.*/,"",$2);print}' |awk '$2=NR $2 $4' > geno_pred/presence_absence_combined/plink/annotated_output_MAF250.bim
cp geno_pred/plink/annotated_output_MAF250.bed geno_pred/presence_absence_combined/plink/annotated_output_MAF250.bed
cp geno_pred/plink/annotated_output_MAF250.fam geno_pred/presence_absence_combined/plink/annotated_output_MAF250.fam


plink --bfile geno_pred/presence_absence_combined/plink/presence_absence_all --bmerge geno_pred/presence_absence_combined/plink/annotated_output_MAF250 --allow-extra-chr --allow-no-sex --make-bed --out geno_pred/presence_absence_combined/plink/annotated_output_MAF250_presence_absence_all


cp geno_pred/presence_absence/Snakefile geno_pred/presence_absence_combined/