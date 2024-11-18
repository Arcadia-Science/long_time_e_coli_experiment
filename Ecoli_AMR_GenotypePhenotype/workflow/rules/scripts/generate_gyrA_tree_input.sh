#!/bin/bash

#################################
input_vcf=$1 #input vcf of all SNPs
gyrA_vcf=$2 #output vcf of all SNPs in gyrA gene
gyrA_vcf_edited=$3 #output vcf of all SNPs in gyrA tree edited so missing genotype is coded as reference
#gyrA_alignment=$4

#select gyrA gene sites only from vcf
cat <(zgrep -e '#' $input_vcf) <(zgrep -w -e 'LMHECDEF_04343' $input_vcf | zgrep -Ev '#') > $gyrA_vcf

#filter out sites with no polymorphism and edit missing genotypes ot reference
bcftools filter -i 'INFO/AC > 0' $gyrA_vcf |\
sed 's/\.\:\./0\:255\,0/g' > $gyrA_vcf_edited

#convert gyrA edited vcf to alignment file for phylogeny construction
 ../vcf2phylip/vcf2phylip.py --output-folder ../tree/alignment/ -i $gyrA_vcf_edited

