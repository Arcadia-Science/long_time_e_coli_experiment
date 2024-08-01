library(snpR)
library(dplyr)
library(ggplot2)


metadata <- read.csv("strains_metadata.csv")
pangenome_samps <- read.csv("pangenome/pangenome_genomes_SRA_GCA.csv") %>% select(SRA) %>% rename(SRA.Accession = SRA)

metadata <- metadata %>% select(SRA.Accession, Isolation.Country, Collection.Year, Season) %>%
        rename(isolation_country = Isolation.Country, collection_year = Collection.Year) %>%
        mutate(isolation_country = gsub(' ','',isolation_country)) %>%
        bind_rows(.,pangenome_samps) %>%
        mutate(sampID = gsub('\\,.*','',SRA.Accession)) %>% select(-SRA.Accession)

dat1 <- import.snpR.data('vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_subsample_head.vcf')
metacomp <- sample.meta(dat1)
missing_meta = subset(sample.meta(dat1) , !(sampID %in% metadata$sampID)) %>% select(sampID)

metadata2 <- metadata %>%
        bind_rows(.,missing_meta)





dat <- import.snpR.data('vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.vcf', sample.meta = metadata2)
#dat <- import.snpR.data('vcf_files/annotated_output_biallelic_synonymous_AC4_refedit_goodcontigs.vcf', sample.meta = metadata2)


p <- plot_structure(dat, facet = 'isolation_country', k = 10, clumpp = FALSE, method = "snapclust")
ggsave('snapclust_country_k10.png',plot = p[[1]], width = 20)

p2 <- plot_structure(dat, facet = 'collection_year', k = 5, clumpp = FALSE, method = "snapclust")
ggsave('snapclust_year.png',plot = p2[[1]], width = 20)







df_clusters <- p[[3]]

df_clusters_k5 <- df_clusters %>% filter(K == 'K = 5') %>% filter( Percentage > 0.9) %>% left_join(metadata)

tabsum <- df_clusters_k5 %>% group_by(Isolation.Country) %>%
mutate(country_n = n()) %>%
group_by(Isolation.Country, Cluster,country_n) %>%
 summarize(count_k = n()) %>% mutate(country_freq = count_k/country_n)












p <- plot_clusters(dat, plot_type = 'tsne')
a <- p$data$pca

tsne <- p$data$
x<- calc_pairwise_ld(dat)
plotld <- plot_pairwise_ld_heatmap(dat)





metadata <- read.csv("strains_metadata.csv")


df <- metadata %>% select(SRA.Accession, Isolation.Country, Size, Season, Sequencing.Platform) %>%
        rename(sampID = SRA.Accession) %>%
        right_join(.,a, by = 'sampID')


pl_pca <- ggplot2::ggplot(df, aes(PC1, PC2, colour = Isolation.Country)) + geom_point()


ggsave('PCA_snpr.png', pl_pca)
