library(data.table)
library(dplyr)
library(tidyr)


#head(geno2[,1:7])


#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
#convert genotype values for rrBLUP to conform to -1 = 0/0
#convert_gt_value <- function(x){ gsub("0", "-1", x)}
#convert_alt_value <- function(x){ gsub("2|3", "1", x)}


#input data
pheno <- fread('strains_metadata_phenotypes_full.txt')

geno <- fread('geno_pred/annotated_output_all_target_genes_tophits_geno.txt')




######################################################
#prep pheno

pheno2 <- pheno %>% select(sample.id, ciprofloxacin, trimethoprim.sulfamethoxazole)



######################################################
#prep geno


#clean column names from bcftools output and create relevant columns
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP) #%>%
                #mutate(across(!snp_id, convert_gt_value))%>%
                #mutate(across(!snp_id, convert_alt_value))


#tranpose to proper orientation and remove duplicate sites
geno_rr_2 <- geno_rr%>%
  group_by(snp_id) %>%
  mutate(count_site = n()) %>% filter(count_site < 2) %>% select(-count_site) %>%
    pivot_longer(cols=c(-snp_id),names_to="sample.id")%>%
    pivot_wider(names_from=c(snp_id))

#head(geno_rr_2[,1:7])

######################################################
#head(pheno_geno[,1:7])
pheno_geno <- inner_join(pheno2, geno_rr_2) #%>% filter(!is.na(ciprofloxacin))


table(pheno_geno$ciprofloxacin, pheno_geno$LMHECDEF_04343_248)
#top hit trim sulfa
table(pheno_geno$trimethoprim.sulfamethoxazole, pheno_geno$`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22`)
#second hit trim sulfa
table(pheno_geno$trimethoprim.sulfamethoxazole, pheno_geno$DHJNCGMO_04404_672)





aminoglycosides <- pheno %>% select(contains('icin'), contains('ycin'))

aminoglycosides <- aminoglycosides %>%
            mutate(across(everything(),~ gsub("Resistant",1, .))) %>%
            mutate(across(everything(),~ gsub("Intermediate",0, .))) %>%
            mutate(across(everything(),~ gsub("Susceptible",-1, .))) %>%
            mutate(across(everything(),~ as.numeric(.)) )

table(aminoglycosides$azithromycin)

corr_pheno <- cor(aminoglycosides, use = 'complete.obs')

cor(aminoglycosides$gentamicin, aminoglycosides$azithromycin, use = 'complete.obs')

cor(aminoglycosides$gentamicin, aminoglycosides$tobramycin, use = 'complete.obs')

table(aminoglycosides$gentamicin, aminoglycosides$tobramycin)

table(aminoglycosides$gentamicin, aminoglycosides$streptomycin)
