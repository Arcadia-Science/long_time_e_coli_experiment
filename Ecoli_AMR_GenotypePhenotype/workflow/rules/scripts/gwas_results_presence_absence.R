library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


#snakemake args
input_gwas_results <- snakemake@input[['gwas_results_combined']]
input_gene_name <- snakemake@input[['gene_names']]

output_tophits_fig <- snakemake@output[['gwas_hits_fig']]
gwas_top10_table <- snakemake@output[['gwas_top10_table']]
gwas_top10_plink_sites <- snakemake@output[['gwas_top10_plink_sites']]
gwas_top10_sites <- snakemake@output[['gwas_top10_sites']]


#read in gwas output and gene names scraped from pangenome fasta headers
df <- fread(input_gwas_results)
gene_names <- fread(input_gene_name, fill = TRUE, header = FALSE) %>% rename(chr = V1, gene = V2)


#merge genomic prediction results to gene names
df  <- df %>% arrange(-gamma) %>% rename_at(ncol(df), ~"phenotype") %>% left_join(., gene_names)


#calculate absolute effect size and sort markers from largest effect to smallest, then label presence/absence loci vs. SNPs based on bogus position coordinate 696969
#effect size estimate is alpha + beta*gamma as per GEMMA recommendation
df <- df %>% mutate(eff =  abs(beta*gamma + alpha))  %>% arrange(-eff) %>%
        mutate(type_marker = ifelse(ps == 696969, "pres_abs", "snp"))


#filter down to top 10 markers (by absolute effect size) per phenotype
top_hits <- df %>% group_by(phenotype) %>%
    mutate(rank = rank(-eff)) %>%
    group_by(phenotype, type_marker) %>%
    mutate(rank_within = rank(-eff)) %>%
        filter(rank < 11) %>% arrange(-eff) %>% data.frame(.)


#print table of annotated genes in top 10 markers
table(top_hits$gene)


#plot marker effect size decay curves for top 10 marker by phenotype
pl1 <- ggplot(data = top_hits, aes(x = rank, y = eff)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        geom_line() + geom_point(aes(x = rank, y = eff, colour = type_marker)) + facet_wrap(~phenotype) +
        xlab('Marker rank') + ylab('Marker effect size') + labs(color = 'Marker type') +
        scale_x_continuous(breaks = function(rank) unique(floor(pretty(seq(min(rank), (max(rank) + 1) * 1.1))))) +
        scale_color_manual(labels = c("Presence/Absence", "SNP"), values = c("#5088C5", "#F28360"))

ggsave(output_tophits_fig, pl1, width = 9, height = 3)


#output table of top 10 markers (by effect size) for each phenotype
write.table(top_hits, gwas_top10_table, row.names = FALSE, quote = FALSE, sep = '\t')


#output sites for plink subsetting of top 10 markers
plink_sites <- top_hits %>% select(chr, ps, phenotype) %>%
        mutate(ps2 = ps, rangeid = paste0(phenotype, ps)) %>%
        select(chr, ps, ps2, rangeid)

write.table(plink_sites, gwas_top10_plink_sites, row.names = FALSE, quote = FALSE, sep = '\t')


#output sites for bcftools subsetting of top 10 markers
bcftools_sites <- top_hits %>% select(chr, ps)
write.table(bcftools_sites, gwas_top10_sites, row.names = FALSE, quote = FALSE, sep = '\t', col.names = FALSE)
