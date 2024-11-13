#!/bin/bash

#################################
input_vcf=$1
gyrA_vcf=$2
gyrA_vcf_edited=$3
#gyrA_alignment=$4

#gyrA gene
cat <(zgrep -e '#' $input_vcf) <(zgrep -w -e 'LMHECDEF_04343' $input_vcf | zgrep -Ev '#') > $gyrA_vcf

bcftools filter -i 'INFO/AC > 0' $gyrA_vcf |\
sed 's/\.\:\./0\:255\,0/g' > $gyrA_vcf_edited

 ../vcf2phylip/vcf2phylip.py --output-folder ../tree/alignment/ -i $gyrA_vcf_edited

#iqtree -s tree/alignment_gyra/annotated_output_gyrA_refedit.min4.phy  -m GTR+ASC
