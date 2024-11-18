library(data.table)
library(dplyr)
library(tidyr)

#snakemake input
input_file <- snakemake@input[[1]]
metadata_pheno_file <- snakemake@input[[2]]
pheno_select <- snakemake@params[[1]]

#fam <- fread('plink/annotated_output_MAF250.fam') %>% rename(sample.id = V1 )
#df <- fread('../strains_metadata_phenotypes_full.txt')


print(pheno_select) #print focal phenotype
df <- fread(metadata_pheno_file) #read in metadata + phenotypes
fam <- fread(input_file) %>% rename(sample.id = V1 ) #read in plink fam file with no phenotypes



phenos <- df %>% select(one_of('sample.id' ,pheno_select)) #select focal phenotype and sample ID from metadata
#recode resistance phenotype to binary state, combining resistant/intermediate together
phenos2 <- phenos %>%
            mutate(across(everything(),~ gsub("Resistant",1, .))) %>%
            mutate(across(everything(),~ gsub("Intermediate",1, .))) %>%
            mutate(across(everything(),~ gsub("Susceptible",0, .))) %>%
            mutate(across(!sample.id,~ as.numeric(.)) )
#join recoded phenotype of interest to fam file
fam_full <- fam %>% select(-V6) %>% left_join(.,phenos2)

write.table(fam_full, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
