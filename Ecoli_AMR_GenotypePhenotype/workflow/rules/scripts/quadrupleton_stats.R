library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#df <- fread('../vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt')
input_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]

#read in file of genotypes for all quadrupleton synonymous sites
df <- fread(input_file)

#funs to clean bcftools output
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
convert_numeric <- function(x){ as.numeric(x)}

#take genotypes and clean column names (to obtain sample names), convert genotypes to numeric
df2 <- df %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP, -V6) %>%
                mutate(across(!snp_id, convert_numeric))

#calculate sums of qudrupletons for all samples (Columns)
df3 <- df2 %>% select(-snp_id) %>% colSums(., na.rm=T)

#pull out outliers and write their sample names to file
outlier_samples <- df3 %>% as.data.frame(.) %>% filter(`.` > 500) %>% tibble:::rownames_to_column(.) %>% select(rowname)

write.table(outlier_samples, out_file, quote = F, col.names = F, row.names = F, sep = '\t' )

#write.table(outlier_samples, 'pangenome/outlier_samples.txt', quote = F, col.names = F, row.names = F, sep = '\t' )