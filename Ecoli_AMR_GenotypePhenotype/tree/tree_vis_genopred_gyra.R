library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)

tree_test<- read.newick ("tree/alignment_gyra/annotated_output_gyrA_refedit.min4.phy.treefile")

metadata <- fread("strains_metadata_phenotypes_full.txt")
geno <- fread('geno_pred/annotated_output_all_target_genes_geno.txt')

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


#simplyify MLST groups to aid in plotting
metadata <- metadata %>% group_by(Other.Typing) %>%
                mutate(count_type = n()) %>%
                mutate(mlst = ifelse(count_type < 100, 'NA', Other.Typing)) %>% ungroup()



pl <- ggtree(tree_test,layout='circular') %<+% metadata
ggsave('figs/test_tree_iqtree_gyra_gamma.png', plot = pl)


pl4 <- pl + geom_tippoint(aes(colour = mlst),size = 0.001) +  theme(legend.position = "none")
ggsave('figs/test_tree_iqtree_gyra_mlst.png', plot = pl4)
#pl4_2 <- pl + geom_tippoint(aes(colour = mlst),size = 0.001)
#ggsave('figs/test_tree_mlst_legend.png', plot = pl4_2)



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
ggsave('figs/test_tree_iqtree_gyra_antibiotic.png', plot = pl5)


##############################
#markers

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                filter(CHROM == 'LMHECDEF_04343') %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248")

geno_rr_subset <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(LMHECDEF_04343_259 = as.factor(LMHECDEF_04343_259),
           LMHECDEF_04343_248 = as.factor(LMHECDEF_04343_248))


metadata_subset <- metadata %>% select( taxa, ciprofloxacin) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin))



geno_rr_subset2 <- geno_rr_subset %>% left_join(.,metadata_subset) %>% select(-taxa) %>% as.data.frame()

rownames(geno_rr_subset2) <- geno_rr_subset$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, geno_rr_subset2,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none")
ggsave('figs/test_tree_iqtree_gyra_marker_gamma.png', plot = pl7)



#plot effects of markers
plot_markers <- geno_rr_subset2 %>% mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                        mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(LMHECDEF_04343_248,LMHECDEF_04343_259))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw()

ggsave('figs/histogram_ciprofloxacin_markers.png', plot = pl_effect, width = 5, height = 5)
