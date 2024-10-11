library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)

df <- fread('geno_pred/presence_absence_combined/output/presence_absence_ALL_bslmm_blup.param.txt')
df  <- df %>% arrange(-gamma) %>% rename_at(ncol(df), ~"phenotype" ) %>% left_join(.,gene_names)
#df <- fread('output/gwas_MAF500_bslmm_probit.param.txt') %>% arrange(-gamma) %>% left_join(.,gene_names)
#head(df)

#df$eff <- df$beta*df$gamma
df <- df %>% mutate(eff =  abs(beta*gamma + alpha))  %>% arrange(-eff) %>% mutate (type_marker = ifelse(ps == 696969, "pres_abs","snp"))
head(df)

#head(df[df$chr == 'LMHECDEF_04343',])

top_hits <- df %>% group_by(phenotype) %>%
    mutate(rank = rank(-eff))%>%
    group_by(phenotype, type_marker) %>%
    mutate(rank_within = rank(-eff)) %>%
        filter(rank <11) %>% arrange(-eff) %>% data.frame(.)



table(top_hits$gene)



pl1 <- ggplot(data = top_hits, aes(x = rank, y = eff)) +
        #geom_line(aes(y = cumsum(eff))) + facet_wrap(~phenotype)
        geom_line() + geom_point( aes(x = rank, y = eff, colour = type_marker)) +  facet_wrap(~phenotype)


ggsave('figs/gwas_presence_absence_complete_bslmm_probit_top10_marker_effects.png', pl1, width = 10, height = 5)


#write.table(top_hits, 'geno_pred/gwas_presence_absence_complete_bslmm_probit_top10.txt', row.names = F, quote = F, sep = '\t')

plink_sites <- top_hits %>% select(chr, ps, phenotype)%>%
        mutate(ps2 = ps, rangeid = paste0(phenotype, ps)) %>%
        select(chr, ps, ps2, rangeid)
#write.table(plink_sites, 'gwas_presence_absence_complete_bslmm_probit_top10_ranges.txt', row.names = F, quote = F, sep = '\t')

bcftools_sites <- top_hits %>% select(chr, ps)
#write.table(bcftools_sites, 'gwas_presence_absence_complete_bslmm_probit_top10_sites.txt', row.names = F, quote = F, sep = '\t', col.names = F)
