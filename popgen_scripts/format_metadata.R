library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)



metadata <- read.csv("strains_metadata.csv")
pheno <- read.csv('phenotype_matrix/phenotype_matrix_08302024.csv')
pheno_name <- read.csv('phenotype_matrix/SRA_to_genome_name.csv')
pangenome_samps <- read.csv("pangenome/pangenome_genomes_SRA_GCA.csv") %>% select(SRA) %>%
        rename(sample.id = SRA) %>% mutate(ref = 'ref')


#right join phenotypes to file mapping names as names file has all strains and pheno is missing two
#then attach pangenome samples to bottom for full list of strains
pheno_dat <- pheno %>% right_join(.,pheno_name, by = 'Genome.Name') %>%
        rename(sample.id = SRA.Accession) %>% mutate(ref = 'non-ref') %>%
        bind_rows(.,pangenome_samps) %>%
        select(-X.x,-X.y)

#join phenotype data to metadata
full_pheno_metadata <- metadata %>%
         mutate(sample.id = gsub('\\,.*','',SRA.Accession)) %>%
        right_join(.,pheno_dat, by = 'sample.id')



write.table(full_pheno_metadata, 'strains_metadata_phenotypes_full.txt', sep = '\t', row.names = F, quote = F)
