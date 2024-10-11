library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)



tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")
#tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.timetree.nwk")


metadata <- fread("strains_metadata_phenotypes_full.txt")

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


#simplyify MLST groups to aid in plotting
metadata <- metadata %>% group_by(Other.Typing) %>%
                mutate(count_type = n()) %>%
                mutate(mlst = ifelse(count_type < 100, '', Other.Typing))


#cluster_data <- fread("tables/snapclust_synonymous_goodcontigs_k2_k10_results.txt")


pl <- ggtree(tree_test,layout='circular') %<+% metadata
ggsave('figs/test_tree_iqtree_remoutliers_mac10.png', plot = pl)

pl2 <- ggtree(tree_test,layout='circular') + geom_text(aes(label=node), hjust=-.3,size=0.75)
ggsave('figs/test_tree_iqtree_remoutliers_mac10_nodelabels.png', plot = pl2)


#pl2 <- pl + geom_tippoint(aes(colour = Isolation.Country),size = 0.0001)
#ggsave('figs/test_tree_country.png', plot = pl2)

#pl3 <- pl + geom_tippoint(aes(colour = as.factor(collection_year)),size = 0.001)
#ggsave('figs/test_tree_year.png', plot = pl3)



###################################################
#annotated trees

#antibiotics
metadata_subset <- metadata %>% ungroup(.) %>% select( ciprofloxacin, ampicillin,gentamicin, cefotaxime) %>%
        mutate_if(is.character, as.factor)


metadata_subset <- as.data.frame(metadata_subset)

rownames(metadata_subset) <- metadata$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl5 <- gheatmap(pl4, metadata_subset,width=0.5,color = NULL,colnames_angle=45, font.size=3) +
    scale_fill_viridis_d(option="D")
ggsave('figs/test_tree_iqtree_remoutliers_mac10_antibiotic.png', plot = pl5)


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
metadata_subset <- metadata %>% ungroup(.) %>% select( mlst) %>%
        mutate(mlst = as.factor(mlst))

metadata_subset <- as.data.frame(metadata_subset)

rownames(metadata_subset) <- metadata$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, metadata_subset,width=0.2,color = NULL,font.size=0) + theme(legend.position = "right")
ggsave('figs/test_tree_iqtree_remoutliers_mac10_mlst.png', plot = pl7)


#mlst plus clade labels
pl7_2 <- gheatmap(pl4, metadata_subset,width=0.2,color = NULL,font.size=0) + theme(legend.position = "right")+
        geom_cladelabel(node=7113, label="phylogroup B2",color="purple",offset=.2,,offset.text=.1) +
        geom_cladelabel(node=11593, label="phylogroup D",color="olivedrab",offset=.33,offset.text=.1) +
        #geom_cladelabel(node=12489, label="phylogroup A",color="blue",offset=.2) +
        geom_cladelabel(node=11297, label="phylogroup F",color="hotpink2",offset=.3,offset.text=.1, hjust=0.8)

ggsave('figs/test_tree_iqtree_remoutliers_mac10_mlst_cladelab.png', plot = pl7_2)
