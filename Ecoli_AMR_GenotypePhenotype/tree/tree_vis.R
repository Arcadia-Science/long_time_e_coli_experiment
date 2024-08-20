library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)



tree_test <- read.newick("tree/alignment/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample.min4.phy.treefile")

#tree_test<- read.newick ("tree/alignment_gyra/annotated_output_gyrA_refedit.min4.phy.treefile")

metadata <- fread("strains_metadata_phenotypes_full.txt")

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


#simplyify MLST groups to aid in plotting
metadata <- metadata %>% group_by(Other.Typing) %>%
                mutate(count_type = n()) %>%
                mutate(mlst = ifelse(count_type < 100, 'NA', Other.Typing))


#cluster_data <- fread("tables/snapclust_synonymous_goodcontigs_k2_k10_results.txt")


pl <- ggtree(tree_test,layout='circular') %<+% metadata
ggsave('figs/test_tree_iqtree_gyra.png', plot = pl)

pl2 <- pl + geom_tippoint(aes(colour = Isolation.Country),size = 0.0001)
ggsave('figs/test_tree_country.png', plot = pl2)

pl3 <- pl + geom_tippoint(aes(colour = as.factor(collection_year)),size = 0.001)
ggsave('figs/test_tree_year.png', plot = pl3)


pl4 <- pl + geom_tippoint(aes(colour = mlst),size = 0.001) +  theme(legend.position = "none")
ggsave('figs/test_tree_iqtree_gyra_mlst.png', plot = pl4)
pl4_2 <- pl + geom_tippoint(aes(colour = mlst),size = 0.001)
ggsave('figs/test_tree_mlst_legend.png', plot = pl4_2)



###################################################
#annotated trees

#antibiotics
metadata_subset <- metadata %>% select( ciprofloxacin, ampicillin) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin),
                ampicillin = as.factor(ampicillin))

metadata_subset <- as.data.frame(metadata_subset)

rownames(metadata_subset) <- metadata$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl5 <- gheatmap(pl4, metadata_subset,width=0.2,color = NULL,colnames_angle=45) +
    scale_fill_viridis_d(option="D")
ggsave('figs/test_tree_cat.png', plot = pl5)


##############################
#clustering groups
cluster2 <- cluster_data %>% select(Cluster) %>% mutate(Cluster = as.factor(Cluster))
cluster2 <- as.data.frame(cluster2)

rownames(cluster2) <- cluster_data$sampID
pl4 <- ggtree(tree_test,layout='circular')

pl6 <- gheatmap(pl4, cluster2,width=0.2,color = NULL,colnames_angle=45) +
    scale_fill_viridis_d(option="D")

ggsave('figs/test_tree_clust.png', plot = pl6)


##############################
#phylogroup
metadata_subset <- metadata %>% select( mlst) %>%
        mutate(mlst = as.factor(mlst))

metadata_subset <- as.data.frame(metadata_subset)

rownames(metadata_subset) <- metadata$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, metadata_subset,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none")
ggsave('figs/test_tree_mlst_cat.png', plot = pl7)
