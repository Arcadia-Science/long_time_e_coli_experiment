#!/bin/bash



#linearize reference fasta
sed -e 's/\(^>.*$\)/#\1#/' pangenome/whole_pangenome.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | \
#extract names of contigs in focal vcf file and read those in as a list of strings to grep the linearized fasta
grep -f <(grep -Ev '#' vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_ldpruned_subsample.vcf | awk '{print $1}' ) -w -A 1 | \
#count number of characters to get site number
grep -Ev '>|--' | wc

