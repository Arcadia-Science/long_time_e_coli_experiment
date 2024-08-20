
grep -e '>' pangenome/whole_pangenome.fasta | sed 's/>//g' > pangenome/gene_names.txt