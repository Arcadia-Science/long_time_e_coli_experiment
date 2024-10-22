library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

#df <- fread('geno_pred/plink/TESTBED.epi.cc', header = T) %>% arrange(P)

#gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)


#df2 <- df %>% left_join(.,gene_names, by = c('CHR1' = 'chr')) %>%
#    rename(genechr1 = gene) %>%
#        left_join(.,gene_names, by = c('CHR2' = 'chr')) %>%
#            rename(genechr2 = gene)

#head(df2)



############

df <- fread('geno_pred/plink/presence_absence_complete_probit_tophits.ld', header = T) %>% arrange(-R2)
gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chr = V1, gene = V2)

df2 <- df %>% left_join(.,gene_names, by = c('CHR_A' = 'chr')) %>%
    rename(genechr1 = gene) %>%
        left_join(.,gene_names, by = c('CHR_B' = 'chr')) %>%
            rename(genechr2 = gene)


hits <- fread('geno_pred/gwas_presence_absence_complete_bslmm_probit_top10.txt', header = T)
hits <- hits %>% mutate(snpid = paste0(phenotype, rank)) %>% select(chr, ps, phenotype) %>%
    group_by(chr, ps) %>% summarize(pheno = paste0(phenotype, collapse = ""))

df3 <- df2 %>% left_join(.,hits, by = c('CHR_A' = 'chr', 'BP_A' = 'ps')) %>%
    rename(phenochr1 = pheno) %>%
        left_join(.,hits, by = c('CHR_B' = 'chr', 'BP_B' = 'ps')) %>%
            rename(phenochr2 = pheno) %>%
            mutate(within_pheno = ifelse(phenochr1 == phenochr2, 1, 0)) %>%
            mutate(within_gene = ifelse(CHR_A == CHR_B, 1, 0))



pl1 <- ggplot(data = df3, aes(x = as.factor(within_pheno), y = R2)) +
        geom_boxplot() + facet_wrap(~within_gene)


#ggsave('figs/geno_pred_LD_tophits.png', pl1)


#write.table(df3, 'tables/gwas_presence_absence_complete_bslmm_probit_top10_ld.txt', quote = F, sep = '\t', row.names = F)

#heatmap for ampicillin
top_hits <- fread('geno_pred/gwas_presence_absence_ALL_bslmm_probit_top_hits.txt', header = T) %>%
    select(rs, rank, phenotype) %>% filter(phenotype == 'ampicillin') %>%
    mutate(marker_name = paste("Marker", rank)) %>%
    select(-phenotype,-rank)



dfheat_single <- df3 %>%  filter(phenochr1 == 'ampicillin' & phenochr2 == 'ampicillin') %>% select(SNP_A, SNP_B, R2)
dfheat_mirror <- df3 %>%  filter(phenochr1 == 'ampicillin' & phenochr2 == 'ampicillin') %>% select(SNP_A, SNP_B, R2) %>%
    rename(SNP_B_temp = SNP_A) %>% rename(SNP_A = SNP_B) %>% rename(SNP_B = SNP_B_temp)


dfheat <- rbind(dfheat_single, dfheat_mirror) %>%
    left_join(.,top_hits, by = c('SNP_A' = 'rs')) %>% rename(Marker_A = marker_name) %>%
    left_join(.,top_hits, by = c('SNP_B' = 'rs')) %>% rename(Marker_B = marker_name)



dfheat$Marker_A_ordered <- factor(dfheat$Marker_A, ordered=TRUE, levels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4",
                                                                            "Marker 5", "Marker 6", "Marker 7", "Marker 8",
                                                                            "Marker 9", "Marker 10"))
dfheat$Marker_B_ordered <- factor(dfheat$Marker_B, ordered=TRUE, levels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4",
                                                                            "Marker 5", "Marker 6", "Marker 7", "Marker 8",
                                                                            "Marker 9", "Marker 10"))



plheat <- ggplot(data = dfheat, aes(x=Marker_A_ordered, y=Marker_B_ordered, fill=R2)) +
        geom_tile() + scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        scale_fill_gradientn(limits = c(0,0.75), colours=c("#341E60", "#A96789", "#F5DFB2"))+
        xlab("First marker") + ylab("Second marker")


#ggsave('final_figs/geno_pred_LD_tophits_heatmap.svg', plheat, width = 6.5, height = 5)

##########################
#null distribution

df_null <- fread('geno_pred/presence_absence_combined_ld_marker_null_set.ld')
df_null$within_gene <- ifelse(df_null$CHR_A == df_null$CHR_B, 1, 0)
df_null$CHR_A_presabs <- ifelse(grepl('pres',df_null$SNP_A ), 1, 0)
df_null$CHR_B_presabs <- ifelse(grepl('pres',df_null$SNP_B ), 1, 0)


df_null_interchr_presabs <- df_null %>% filter(within_gene == 0 & CHR_A_presabs == 1 & CHR_B_presabs == 1)
df_null_interchr_snps <- df_null %>% filter(within_gene == 0 & CHR_A_presabs == 0 & CHR_B_presabs == 0)
df_null_interchr_between <- df_null %>% filter(within_gene == 0 & CHR_A_presabs == 0 & CHR_B_presabs == 1 | within_gene == 0 & CHR_A_presabs == 1 & CHR_B_presabs == 0 )
df_null_within_gene <- df_null %>% filter(within_gene == 1 & CHR_A_presabs == 0 & CHR_B_presabs == 0)


#calculate p-values for tophits
df_pval <- df3 %>% mutate(CHR_A_presabs = ifelse(grepl('pres',SNP_A ), 1, 0),
                         CHR_B_presabs = ifelse(grepl('pres',SNP_B ), 1, 0)) %>%
                            mutate(type = ifelse(CHR_A_presabs == 0 & CHR_A_presabs == 0, 'snps',
                                            ifelse(CHR_A_presabs == 1 & CHR_A_presabs == 1, 'pres_abs','mixed')))


df_pval <- df_pval %>% mutate(prcnt =
                        ifelse(type == 'snps' & within_gene == 0, ecdf(df_null_interchr_snps$R2)(R2), #inter chr SNPs
                            ifelse(type == 'pres_abs' & within_gene == 0, ecdf(df_null_interchr_presabs$R2)(R2), #inter chr presabs
                                 ifelse(type == 'snps' & within_gene ==1, ecdf(df_null_within_gene$R2)(R2), #inter chr mixed snp/presabs
                                        ifelse(type == 'snps' & within_gene ==1, ecdf(df_null_within_gene$R2)(R2),NA))))) #within gene snp

df_pval$pval <- ifelse(df_pval$prcnt < 0.5, df_pval$prcnt, (1-df_pval$prcnt))
#write.table(df_pval, 'tables/gwas_presence_absence_complete_bslmm_probit_top10_ld.txt', quote = F, sep = '\t', row.names = F)
