library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggcorrplot)



df <- fread('strains_metadata_phenotypes_full.txt')


df2 <- df %>% select(Genome.Name.x, trimethoprim.sulfamethoxazole, Isolation.Country, sample.id)



df2 <- df %>% select(Genome.Name.x, trimethoprim.sulfamethoxazole, Isolation.Country, ciprofloxacin, ampicillin)

df3 <- df2 %>% filter(trimethoprim.sulfamethoxazole== 'Resistant')


df4 <- df2 %>% group_by(trimethoprim.sulfamethoxazole ) %>%
        mutate(country_m = n()) %>%
            group_by(Isolation.Country, trimethoprim.sulfamethoxazole, country_m) %>%
            summarize(coun_pheno = n()) %>% mutate(prop_pheno = coun_pheno/country_m) %>%
            filter(trimethoprim.sulfamethoxazole == 'Resistant') %>% data.frame(.)





df4 <- df2 %>% group_by(ciprofloxacin) %>%
        mutate(country_m = n()) %>%
            group_by(Isolation.Country, ciprofloxacin, country_m) %>%
            summarize(coun_pheno = n()) %>% mutate(prop_pheno = coun_pheno/country_m) %>%
            filter(ciprofloxacin == 'Resistant') %>% data.frame(.)


df4 <- df2 %>% group_by(ampicillin) %>%
        mutate(country_m = n()) %>%
            group_by(Isolation.Country, ampicillin, country_m) %>%
            summarize(coun_pheno = n()) %>% mutate(prop_pheno = coun_pheno/country_m) %>%
            filter(ampicillin == 'Resistant') %>% data.frame(.)





################################################################################################################
################################################################################################################
################################################################################################################
metadata <- fread("strains_metadata_phenotypes_full.txt")
geno <- fread('geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt')

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


#simplyify MLST groups to aid in plotting
metadata <- metadata %>% group_by(Other.Typing) %>%
                mutate(count_type = n()) %>%
                mutate(mlst = ifelse(count_type < 100, 'NA', Other.Typing)) %>% ungroup()



#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                #filter(CHROM == 'LMHECDEF_04343') %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22" | snp_id == "DHJNCGMO_04404_672")

geno_rr_subset <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22` = as.factor(`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22`),
           DHJNCGMO_04404_672 = as.factor(DHJNCGMO_04404_672))


metadata_subset <- metadata %>% select( taxa, trimethoprim.sulfamethoxazole, Isolation.Country) %>%
        mutate(trimethoprim.sulfamethoxazole = as.factor(trimethoprim.sulfamethoxazole))



geno_rr_subset2 <- geno_rr_subset %>% left_join(.,metadata_subset) %>% select(-taxa) %>% as.data.frame() %>%
    mutate(norway = ifelse(Isolation.Country == 'Norway',1,0))


table(geno_rr_subset2$trimethoprim.sulfamethoxazole,geno_rr_subset2$`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22`, geno_rr_subset2$norway)
