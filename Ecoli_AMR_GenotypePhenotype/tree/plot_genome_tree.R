library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(gridExtra)





#tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")
#metadata <- fread("strains_metadata_phenotypes_full.txt")


###########################
tree_file <- snakemake@input[['tree_genome']]
metadata_file <- snakemake@input[['metadata_formatted']]

tree_fig_file <- snakemake@output[['tree_genome_plot']]


###########################
#read input file (tree and metadata)

tree_test <- read.newick(tree_file)
metadata <- fread(metadata)



###########################
#format metadata to ignore small MLST groups (less than 100 samples) that clutter figures
metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )%>%
                 group_by(Other.Typing) %>%
                        mutate(count_type = n()) %>%
                        mutate(mlst = ifelse(count_type < 100, '', Other.Typing))


###########################
#plot annotated trees

###########
#phylogroup
metadata_subset_phylo <- metadata %>% ungroup(.) %>% select( mlst) %>%
        mutate(mlst = as.factor(mlst)) %>% data.frame(.)

rownames(metadata_subset_phylo) <- metadata$taxa

#mlst plus clade labels
#mlst plus clade labels
plot_phylogroups <- gheatmap(ggtree(tree_test,layout='circular'), metadata_subset_phylo,width=0.2,color = NULL,font.size=0) +
        theme(legend.position = "right") +
        geom_cladelabel(node=7113, label="phylogroup B2",color="purple",offset=.2,,offset.text=.1) +
        geom_cladelabel(node=11593, label="phylogroup D",color="olivedrab",offset=.33,offset.text=.1) +
        geom_cladelabel(node=11297, label="phylogroup F",color="hotpink2",offset=.3,offset.text=.1, hjust=0.8)





plot_phylogroups_test <- gheatmap(ggtree(tree_test,layout='circular'), metadata_subset_phylo,width=0.2,color = NULL,font.size=0) +
        theme(legend.position = "right", legend.text=element_text(size=1), legend.key.size = unit(0.05, 'cm'))+
        geom_cladelabel(node=7113, label="phylogroup B2",color="purple",offset=.2,,offset.text=.1) +
        geom_cladelabel(node=11593, label="phylogroup D",color="olivedrab",offset=.33,offset.text=.1) +
        geom_cladelabel(node=11297, label="phylogroup F",color="hotpink2",offset=.3,offset.text=.1, hjust=0.8)

#ggsave('figs/phyloforAudrey.pdf', plot_phylogroups)
#ggsave('figs/phyloforAudrey.pdf', plot_phylogroups_test, width = 430, height = 430, units = 'px')
###########
#antibiotics
metadata_subset <- metadata %>% ungroup(.) %>% select( ciprofloxacin, ampicillin, trimethoprim.sulfamethoxazole) %>%
        mutate_if(is.character, as.factor) %>% data.frame(.)

rownames(metadata_subset) <- metadata$taxa

#focal AMR phenotype states
plot_phenotypes <- gheatmap(ggtree(tree_test,layout='circular'), metadata_subset,width=0.5,color = NULL,colnames_angle=45, font.size=3) +
    scale_fill_viridis_d()

###########
#combine plots
#combine plots into one plot
plot_trees_all <- grid.arrange(plot_phylogroups, plot_phenotypes, ncol=2, nrow =1)

#output figure
ggsave(tree_fig_file, plot_trees_all, width = 16)
