library(SNPRelate)
library(ggplot2)
library(dplyr)


vcf.fn <- "vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_subsample.vcf"

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")


genofile <- snpgdsOpen("test.gds", readonly = F)
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)


metadata <- read.csv("strains_metadata.csv")
pheno <- read.csv('phenotype_matrix/phenotype_matrix_08302024.csv')
pheno_name <- read.csv('phenotype_matrix/SRA_to_genome_name.csv')

pheno_dat <- pheno %>% left_join(.,pheno_name, by = 'Genome.Name') %>%
        select(SRA.Accession, ciprofloxacin ) %>%
        rename(sample.id = SRA.Accession)

df <- metadata %>% select(SRA.Accession, Isolation.Country, Size, Season, Sequencing.Platform) %>%
        rename(sample.id = SRA.Accession) %>%
        right_join(.,tab, by = 'sample.id') %>%
        left_join(.,pheno_dat, by = 'sample.id')


pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()

#pl <- ggplot(df, aes(x=EV1, y=EV2, colour = ciprofloxacin)) + geom_point()



df2 <- df  %>% filter(Isolation.Country != 'Norway' & Isolation.Country != 'United Kingdom'& Isolation.Country != 'England')
pl <- ggplot(df2, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()+ stat_ellipse()

pl <- ggplot(df2, aes(x=EV1, y=EV2, colour = ciprofloxacin)) + geom_point()+ stat_ellipse()


ggsave('PCA.pdf', pl)
