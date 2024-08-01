library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)




tree_test <- read.newick("tree/alignment/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.min4.phy.treefile")


metadata <- read.csv("strains_metadata.csv")
pangenome_samps <- read.csv("pangenome/pangenome_genomes_SRA_GCA.csv") %>% select(SRA) %>% rename(SRA.Accession = SRA)
metadata <- metadata %>% select(SRA.Accession, Isolation.Country, Collection.Year, Season) %>%
        rename(isolation_country = Isolation.Country, collection_year = Collection.Year) %>%
        mutate(isolation_country = gsub(' ','',isolation_country)) %>%
        bind_rows(.,pangenome_samps) %>%
        mutate(taxa = gsub('\\,.*','',SRA.Accession)) %>% select(-SRA.Accession) %>%
        relocate(taxa, .before = isolation_country )



pl <- ggtree(tree_test,layout='circular') %<+% metadata
ggsave('test_tree.png', plot = pl)

pl2 <- pl + geom_tippoint(aes(colour = isolation_country),size = 0.0001)
ggsave('test_tree_country.png', plot = pl2)

pl3 <- pl + geom_tippoint(aes(colour = as.factor(collection_year)),size = 0.001)
ggsave('test_tree_year.png', plot = pl3)
