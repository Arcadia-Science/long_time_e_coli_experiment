library(SNPRelate)
library(ggplot2)
library(dplyr)
library(data.table)
library(gridExtra)

vcf.fn <- "vcf_files/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.vcf"

snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")


genofile <- snpgdsOpen("test.gds", readonly = F)
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)


metadata <- fread("strains_metadata_phenotypes_full.txt")


df <- metadata %>% select(sample.id, Isolation.Country, Size, Season, Sequencing.Platform, Collection.Year,ref, Other.Typing) %>%
        right_join(.,tab, by = 'sample.id') %>%
                group_by(Other.Typing) %>% mutate(ntype = n()) %>%
                mutate(mlst = ifelse(ntype > 49, Other.Typing, NA) ) %>% ungroup(.)


pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Isolation.Country)) + geom_point()
#pl <- ggplot(df, aes(x=EV1, y=EV2, colour = Collection.Year)) + geom_point()

ggsave('figs/PCA_goodcontigs_remoutliers_country.png', pl)

##############
#remove larger country clusters to aid in datavis
#df2 <- df  %>% filter(Isolation.Country != 'Norway' & Isolation.Country != 'United Kingdom'& Isolation.Country != 'England')

#other viz
pl_ref <- ggplot(df, aes(x=EV1, y=EV2, colour = ref)) + geom_point()
pl_mlst <-ggplot(df, aes(x=EV1, y=EV2, colour = mlst)) + geom_point()+ stat_ellipse()
pl_year <-ggplot(df, aes(x=EV1, y=EV2, colour = as.character(Collection.Year))) + geom_point() + labs(colour='Year')


pl2 <- grid.arrange(pl_ref, pl_year, pl_mlst, ncol=3, nrow =1)

ggsave('figs/PCA_goodcontigs_remoutliers_combined_figs.png', pl2, width = 21)

#ggsave('figs/PCA_goodcontigs_remoutliers_mlst.png', pl_mlst, width = 13, height = 13)
