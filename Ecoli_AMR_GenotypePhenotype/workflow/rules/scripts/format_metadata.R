#!/usr/bin/Rscript


library(data.table)
library(dplyr)
library(tidyr)


##########################
#input file names from snakemake
SRA_GCA_file <- snakemake@input[["SRA_GCA"]]
phenotype_matrix_file <- snakemake@input[["phenotype_matrix"]]
SRA_genome_name_file <- snakemake@input[["SRA_genome_name"]]
strains_metadata_file <- snakemake@input[["strains_metadata"]]
outlier_samples_file <- snakemake@input[["outlier_samples"]]

#output filenames from snakemake
full_pheno_metadata_file <- snakemake@output[["metadata_formatted"]]
samplist_file <- snakemake@output[["sample_list"]]
time_calibration_file <- snakemake@output[["time_calibration_data"]]


##########################
#read in input files

#file of pangenome sample names, read them in, select their SRA ID's and create a new ref column that identifies as pangenome samples
#gsub command to edit cases where multiple SRA id's are included for one sample, just take the first one in the comma separated list
pangenome_samps <- read.csv(SRA_GCA_file) %>% select(SRA) %>%
        mutate(sample.id = gsub("\\,.*", "", SRA)) %>% mutate(ref = "ref") %>% select(-SRA)

# metadata file including AMR phenotypes
pheno <- read.csv(phenotype_matrix_file)
#file including SRA's ID's of all samples except for pangenome samples
pheno_name <- read.csv(SRA_genome_name_file)
#metadata file containing alternate data such as sampling location etc
metadata <- read.csv(strains_metadata_file)
#outliers identified in filtering that need to be removed
outliers <- fread(outlier_samples_file, header = FALSE)


##########################
#manipulate metadata


#right join phenotypes to file mapping names as names file has all strains and pheno is missing two
#then attach pangenome samples to bottom for full list of strains
pheno_dat <- pheno %>% right_join(., pheno_name, by = "Genome.Name") %>%
        mutate(sample.id = gsub("\\,.*", "", SRA.Accession)) %>%
        mutate(ref = "non-ref") %>% select(-SRA.Accession) %>%
        bind_rows(., pangenome_samps) %>%
        select(-X.x, -X.y)

#join phenotype data to metadata
full_pheno_metadata <- metadata %>%
         mutate(sample.id = gsub("\\,.*", "", SRA.Accession)) %>%
        right_join(., pheno_dat, by = "sample.id")

#generate full sample list (excluding outliers)
samplist <- full_pheno_metadata %>% select(sample.id)

#generate file for phylogeny time calibration that includes sample ID's and year of sampling
time_calibration <- full_pheno_metadata %>% select(sample.id, Collection.Year) %>%
        filter(!sample.id %in% outliers$V1) %>% filter(!is.na(Collection.Year))



##########################
#write output files
#full metadata set plus phenotypes
#list of samples
#list of samples and year of sampling if available for downstream time calibration of phylogeny

write.table(full_pheno_metadata, full_pheno_metadata_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(samplist, samplist_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(time_calibration, time_calibration_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
