library(data.table)
library(dplyr)
library(tidyr)

gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)

df <- fread('geno_pred/output/gwas_MAF500_ampicillin_bslmm_blup.param.txt') %>% arrange(-gamma) %>% left_join(.,gene_names)
#df <- fread('output/gwas_MAF500_bslmm_probit.param.txt') %>% arrange(-gamma) %>% left_join(.,gene_names)
#head(df)

#df$eff <- df$beta*df$gamma
df <- df %>% mutate(eff =  abs(beta*gamma + alpha))  %>% arrange(-eff)
head(df)

#head(df[df$chr == 'LMHECDEF_04343',])

hits_sum <- df %>% filter(alpha > 0.001) %>% group_by(chr,gene) %>% summarize(count_hits = n()) %>% arrange(-count_hits)
head(hits_sum)

hits_af <- df %>% filter(alpha > 0.001) %>% arrange(-alpha)
head(hits_af)

df2 <- fread('geno_pred/output/gwas_MAF250_ampicillin_bslmm_blup.param.txt') %>% arrange(-alpha) %>% left_join(.,gene_names)