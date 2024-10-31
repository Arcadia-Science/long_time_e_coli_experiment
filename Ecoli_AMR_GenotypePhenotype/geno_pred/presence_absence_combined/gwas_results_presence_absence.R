library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


#read in gene names scraped from pangenome fasta headers
gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)

#genomic prediction results, merged to gene names
df <- fread('geno_pred/presence_absence_combined/output/presence_absence_ALL_bslmm_blup.param.txt')
df  <- df %>% arrange(-gamma) %>% rename_at(ncol(df), ~"phenotype" ) %>% left_join(.,gene_names)


#calculate absolute effect size and sort markers from largest effect to smallest, then label presence/absence loci vs. SNPs based on bogus position coordinate 696969
df <- df %>% mutate(eff =  abs(beta*gamma + alpha))  %>% arrange(-eff) %>%
        mutate (type_marker = ifelse(ps == 696969, "pres_abs","snp"))


#head(df[df$chr == 'LMHECDEF_04343',])

#filter down to top 10 markers (by effect size) per phenotype
top_hits <- df %>% group_by(phenotype) %>%
    mutate(rank = rank(-eff))%>%
    group_by(phenotype, type_marker) %>%
    mutate(rank_within = rank(-eff)) %>%
        filter(rank <11) %>% arrange(-eff) %>% data.frame(.)



table(top_hits$gene)


#plot marker effect size decay curves for top 10 marker by phenotype
pl1 <- ggplot(data = top_hits, aes(x = rank, y = eff)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        geom_line() + geom_point( aes(x = rank, y = eff, colour = type_marker)) +  facet_wrap(~phenotype) +
        xlab('Marker rank') + ylab('Marker effect size') + labs(color='Marker type') +
        scale_x_continuous(breaks = function(rank) unique(floor(pretty(seq(min(rank), (max(rank) + 1) * 1.1)))))+
        scale_color_manual(labels = c("Presence/Absence", "SNP"), values = c("#5088C5", "#F28360"))

#ggsave('final_figs/Fig4_marker_effects_700.png', pl1, width = 9, height = 3)
#ggsave('final_figs/Fig4_marker_effects_700.svg', pl1, width = 9, height = 3)


#write.table(top_hits, 'geno_pred/gwas_presence_absence_complete_bslmm_probit_top10.txt', row.names = F, quote = F, sep = '\t')

plink_sites <- top_hits %>% select(chr, ps, phenotype)%>%
        mutate(ps2 = ps, rangeid = paste0(phenotype, ps)) %>%
        select(chr, ps, ps2, rangeid)
#write.table(plink_sites, 'gwas_presence_absence_complete_bslmm_probit_top10_ranges.txt', row.names = F, quote = F, sep = '\t')

bcftools_sites <- top_hits %>% select(chr, ps)
#write.table(bcftools_sites, 'gwas_presence_absence_complete_bslmm_probit_top10_sites.txt', row.names = F, quote = F, sep = '\t', col.names = F)
