library(snpR)
library(dplyr)
library(ggplot2)
library(data.table)

metadata <- fread("strains_metadata_phenotypes_full.txt")

metadata2 <- metadata %>% select(sample.id, Isolation.Country, Collection.Year, Season) %>%
        rename(isolation_country = Isolation.Country, collection_year = Collection.Year) %>%
        mutate(isolation_country = gsub(' ','',isolation_country)) %>%
        rename(sampID = sample.id)



dat <- import.snpR.data('vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.vcf', sample.meta = metadata2)
#dat <- import.snpR.data('vcf_files/annotated_output_biallelic_synonymous_AC4_refedit_goodcontigs.vcf', sample.meta = metadata2)


p <- plot_structure(dat, facet = 'isolation_country', k = 2:3, clumpp = FALSE, method = "snapclust")
ggsave('figs/snapclust_country_k10.png',plot = p[[1]], width = 20)

p2 <- plot_structure(dat, facet = 'collection_year', k = 2:3, clumpp = FALSE, method = "snapclust")
ggsave('figs/snapclust_year.png',plot = p2[[1]], width = 20)



p <- plot_structure(dat, facet = 'isolation_country', k = 9:10,
        clumpp = FALSE,
        burnin = 100,
        method = "structure",
        structure_path = '/home/ubuntu/miniconda3/envs/popgenR/bin/structure')


ggsave('figs/structure_country_test_k9_10.png',plot = p[[1]], width = 20)


#########################
#output results

#extract cluster assignment
#df_clusters <- p[[3]]
#filter(K == 'K = 10')

df_clusters_k10 <- p[[3]] %>% filter( Percentage > 0.9)
write.table(df_clusters_k10, 'tables/structure_synonymous_goodcontigs_k9_k10_results.txt', sep = '\t', quote = F, row.names =F)