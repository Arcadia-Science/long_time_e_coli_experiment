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
