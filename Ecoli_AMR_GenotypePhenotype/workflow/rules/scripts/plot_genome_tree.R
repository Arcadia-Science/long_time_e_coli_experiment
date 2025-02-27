library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(gridExtra)





###########################
tree_file <- snakemake@input[['tree_genome']]
metadata_file <- snakemake@input[['metadata_formatted']]

tree_fig_file <- snakemake@output[['tree_genome_plot']]


###########################
#read input file (tree and metadata)

tree_test <- read.newick(tree_file)
metadata <- fread(metadata_file)



###########################
#format metadata to ignore small MLST groups (less than 100 samples) that clutter figures
metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )%>%
                 group_by(Other.Typing) %>%
                        mutate(count_type = n()) %>%
                        mutate(mlst_large = ifelse(count_type < 100, NA, Other.Typing)) %>%
                        mutate(phylogroup = ifelse(grepl("phylo",mlst_large), mlst_large, NA)) %>%
                        mutate(mlst = ifelse(grepl("mlst",mlst_large), mlst_large, NA))


###########################
#plot annotated trees

###########
#phylogroup/MLST
metadata_subset_phylo <- metadata %>% ungroup(.) %>% select( mlst, phylogroup) %>%
        mutate(mlst = as.factor(mlst)) %>%
        mutate(phylogroup = as.factor(phylogroup)) %>%
        data.frame(.)
#add sample names as rownames to allow for plotting of metadata on tree
rownames(metadata_subset_phylo) <- metadata$taxa

#mlst/phylogroup clade labels
clade_colours <- c(
                 'mlst:127' = '#B5BEA4',
                 'mlst:12' = '#F28360',
                 'mlst:131' = '#F7B846',
                 'mlst:69' = '#3B9886',
                 'mlst:73' = '#C85152',
                 'mlst:95' = '#8A99AD',
                 'phylogroup:A' = '#73B5E3',
                 'phylogroup:B2' = '#7A77AB',
                 'phylogroup:D' = '#97CD78',
                 'phylogroup:F' = '#F898AE'
                 )
#plot species tree with clade labels
plot_phylogroups <- gheatmap(ggtree(tree_test,layout='circular', size=0.25), metadata_subset_phylo,width=0.3,color = NULL,font.size=0) +
        theme(legend.position = "right") + guides(fill=guide_legend(title="Phylogroup/MLST")) +
        scale_fill_manual(values= clade_colours, na.translate = F)+
        geom_cladelabel(node=7113, label="phylogroup B2",color="#7A77AB",offset=.23,offset.text=.1) +
        geom_cladelabel(node=11593, label="phylogroup D",color="#97CD78",offset=.36,offset.text=.1) +
        geom_cladelabel(node=11297, label="phylogroup F",color="#F898AE",offset=.33,offset.text=.1, hjust=0.8)




###########
#select specific antibiotic resistance phenotypes to plot from metadata, convert to factor
metadata_subset <- metadata %>% ungroup(.) %>% select(ciprofloxacin, ampicillin, trimethoprim.sulfamethoxazole) %>%
        mutate_if(is.character, as.factor) %>% data.frame(.)

#add sample names as rownames to allow for plotting of metadata on tree
rownames(metadata_subset) <- metadata$taxa

#plot tree with focal resistance phenotypes painted on
plot_phenotypes <- gheatmap(ggtree(tree_test, layout = 'circular', size = 0.25), metadata_subset, width = 0.3, color = NULL, colnames_angle = 45, font.size = 3) +
    viridis::scale_fill_viridis(discrete = TRUE, na.translate = FALSE) + guides(fill = guide_legend(title = "AMR phenotype"))


###########
#combine plots into one plot
plot_trees_all <- grid.arrange(plot_phylogroups, plot_phenotypes, ncol = 2, nrow = 1)

#output figure
ggsave(tree_fig_file, plot_trees_all, width = 13)
