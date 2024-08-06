library(SNPRelate)
library(ggplot2)
library(dplyr)


vcf.fn <- "vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.vcf"

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")


genofile <- snpgdsOpen("test.gds", readonly = F)
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)


metadata <- fread("strains_metadata_phenotypes_full.txt")


df <- metadata %>% select(sample.id, Isolation.Country, Size, Season, Sequencing.Platform, Collection.Year) %>%
        right_join(.,tab, by = 'sample.id') 


pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()
#pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Collection.Year)) + geom_point()


##############
#remove larger country clusters to aid in datavis
df2 <- df  %>% filter(Isolation.Country != 'Norway' & Isolation.Country != 'United Kingdom'& Isolation.Country != 'England')
pl <- ggplot(df2, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()+ stat_ellipse()

pl <- ggplot(df2, aes(x=EV1, y=EV2, colour = ciprofloxacin)) + geom_point()+ stat_ellipse()


ggsave('PCA_goodcontigs_year.png', pl)
