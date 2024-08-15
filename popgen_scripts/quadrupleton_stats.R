library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

df <- fread('vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt')

df_geno <- df %>% select(-(1:6),-ncol(df))
sumqaud <- colSums(df_geno, 1)


count(df_geno,1)

df2 <- apply(X = df_geno, MARGIN = 2, FUN = as.numeric)

df_geno2 <- sapply( df_geno, as.numeric )


df2[] <- lapply(df_geno, function(x) as.numeric(as.character(x)))





#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
#convert genotype values for rrBLUP to conform to -1 = 0/0
convert_gt_value <- function(x){ gsub("0", "-1", x)}

convert_numeric <- function(x){ as.numeric(x)}

df2 <- df %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP, -V6) %>%
                mutate(across(!snp_id, convert_numeric))

#viz results
df3 <- df2 %>% select(-snp_id) %>% colSums(., na.rm=T)
table(df3)

#pull out outliers
df3 %>% as.data.frame(.) %>% filter(`.` > 500)