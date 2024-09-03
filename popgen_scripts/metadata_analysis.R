library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggcorrplot)



df <- fread('strains_metadata_phenotypes_full.txt')


outliers <- df %>% filter(sample.id %in% c('ERR1218638','ERR1218722','ERR4035496','ERR4037257','SRR10271534','SRR10272173','SRR10272272','SRR10272285','SRR10272401','SRR2449168','SRR3588897','SRR3932419','SRR3982215','SRR850803'))
table(outliers$ref)
table(outliers$Isolation.Country)
table(outliers$Other.Typing)
table(outliers$Collection.Year)


phenos <- df %>% select(ciprofloxacin, ampicillin, cefotaxime, ceftazidime, gentamicin, amoxicillin.clavulanic.acid, piperacillin.tazobactam,trimethoprim.sulfamethoxazole)
phenos2 <- phenos %>%
            mutate(across(everything(),~ gsub("Resistant",1, .))) %>%
            mutate(across(everything(),~ gsub("Intermediate",0, .))) %>%
            mutate(across(everything(),~ gsub("Susceptible",-1, .))) %>%
            mutate(across(everything(),~ as.numeric(.)) )


corr_pheno <- cor(phenos2, use = 'complete.obs')

p1 <- ggcorrplot(corr_pheno ) + theme(plot.background = element_rect(fill='white', color='white'))
ggsave('figs/phenotype_correlations.png', p1)




png('figs/phenotype_correlations_heatmap.png', width = 1000, height = 1000,res=499)
#par(mar = c(5, 5, 5, 5))
heatmap(x = corr_pheno, symm = TRUE,cex.axis=0.2, margins = c(5,5), cexRow = 0.33,  cexCol = 0.33)
dev.off()