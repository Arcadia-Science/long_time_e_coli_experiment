library(data.table)
library(dplyr)
library(tidyr)


input_file <- snakemake@input[[1]]
pheno_select <- snakemake@params[[1]]

print(pheno_select)
df <- fread('../strains_metadata_phenotypes_full.txt')
fam <- fread(input_file) %>% rename(sample.id = V1 )

#fam <- fread('plink/annotated_output_MAF250.fam') %>% rename(sample.id = V1 )


phenos <- df %>% select(one_of('sample.id' ,pheno_select))
phenos2 <- phenos %>%
            mutate(across(everything(),~ gsub("Resistant",1, .))) %>%
            mutate(across(everything(),~ gsub("Intermediate",1, .))) %>%
            mutate(across(everything(),~ gsub("Susceptible",0, .))) %>%
            mutate(across(!sample.id,~ as.numeric(.)) )
#join phenotype of interest to fam file
fam_full <- fam %>% select(-V6) %>% left_join(.,phenos2)

write.table(fam_full, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
