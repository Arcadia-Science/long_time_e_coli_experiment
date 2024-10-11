library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#df <- fread('vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt')
input_file <- snakemake@input[[1]]
out_file <- snakemake@out[[1]]


df <- fread(input_file)


#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
convert_numeric <- function(x){ as.numeric(x)}

df2 <- df %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP, -V6) %>%
                mutate(across(!snp_id, convert_numeric))

#viz results
df3 <- df2 %>% select(-snp_id) %>% colSums(., na.rm=T)

#pull out outliers
outlier_samples <- df3 %>% as.data.frame(.) %>% filter(`.` > 500) %>% tibble:::rownames_to_column(.) %>% select(rowname)

write.table(outliers_samples, out_file, quote = F, col.names = F, row.names = T, sep = '\t' )

#write.table(outlier_samples, 'pangenome/outlier_samples.txt', quote = F, col.names = F, row.names = F, sep = '\t' )