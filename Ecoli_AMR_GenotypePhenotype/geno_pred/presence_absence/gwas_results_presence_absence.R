library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)

df <- fread('geno_pred/presence_absence/output/presence_absence_ALL_bslmm_blup.param.txt')
df  <- df %>% arrange(-gamma) %>% rename_at(ncol(df), ~"phenotype" ) %>% left_join(.,gene_names)
#df <- fread('output/gwas_MAF500_bslmm_probit.param.txt') %>% arrange(-gamma) %>% left_join(.,gene_names)
#head(df)

#df$eff <- df$beta*df$gamma
df <- df %>% mutate(eff =  abs(beta*gamma + alpha))  %>% arrange(-eff)
head(df)

#head(df[df$chr == 'LMHECDEF_04343',])

top_hits <- df %>% group_by(phenotype) %>%
    mutate(rank = rank(-eff)) %>% filter(rank <11) %>% arrange(-eff) %>% data.frame(.)



table(top_hits$gene)



pl1 <- ggplot(data = top_hits, aes(x = rank, y = eff)) +
        #geom_line(aes(y = cumsum(eff))) + facet_wrap(~phenotype)
        geom_line() + facet_wrap(~phenotype)


ggsave('figs/gwas_presence_absence_bslmm_blup_top10_marker_effects.png', pl1)


#write.table(top_hits, 'geno_pred/gwas_presence_absence_ALL_bslmm_blup_top_hits.txt', row.names = F, quote = F, sep = '\t')

plink_sites <- top_hits %>% select(chr, ps, phenotype)%>%
        mutate(ps2 = ps, rangeid = paste0(phenotype, ps)) %>%
        select(chr, ps, ps2, rangeid)
#write.table(plink_sites, 'geno_pred/gwas_presence_absence_ALL_bslmm_probit_top_hits_plink_ranges.txt', row.names = F, quote = F, sep = '\t')

bcftools_sites <- top_hits %>% select(chr, ps)
#write.table(bcftools_sites, 'geno_pred/gwas_presence_absence_ALL_bslmm_probit_top_hits_sites.txt', row.names = F, quote = F, sep = '\t', col.names = F)
