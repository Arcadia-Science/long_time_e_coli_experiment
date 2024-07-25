library(SNPRelate)
library(ggplot2)
library(dplyr)


vcf.fn <- "vcf_files/annotated_output_biallelic_synonymous_subsample_refedit.vcf"

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")


genofile <- snpgdsOpen("test.gds", readonly = F)
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)


metadata <- read.csv("strains_metadata.csv")

df <- metadata %>% select(SRA.Accession, Isolation.Country, Size, Season, Sequencing.Platform) %>%
        rename(sample.id = SRA.Accession) %>%
        right_join(.,tab, by = 'sample.id')


pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()

df %>% filter(Isolation.Country != 'Norway' & Isolation.Country != 'United Kingdom'& Isolation.Country != 'England') %>%
ggplot( aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()+ stat_ellipse()


#pl <- ggplot(df, aes(x=EV1, y = Size)) + geom_point()

ggsave('PCA.pdf', last_plot())
