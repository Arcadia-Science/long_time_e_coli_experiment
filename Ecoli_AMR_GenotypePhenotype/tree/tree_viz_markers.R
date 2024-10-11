library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)

tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")


metadata <- fread("strains_metadata_phenotypes_full.txt")
geno <- fread('geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt')

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


#simplyify MLST groups to aid in plotting
metadata <- metadata %>% group_by(Other.Typing) %>%
                mutate(count_type = n()) %>%
                mutate(mlst = ifelse(count_type < 100, 'NA', Other.Typing)) %>% ungroup()



pl <- ggtree(tree_test,layout='circular') %<+% metadata
#ggsave('figs/test_tree_iqtree_gamma.png', plot = pl)




###################################################
#annotated trees

#antibiotics
metadata_subset <- metadata %>% select( ciprofloxacin, ampicillin, trimethoprim.sulfamethoxazole) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin),
                ampicillin = as.factor(ampicillin),
                trimethoprim.sulfamethoxazole = as.factor(trimethoprim.sulfamethoxazole))

metadata_subset <- as.data.frame(metadata_subset)

rownames(metadata_subset) <- metadata$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl5 <- gheatmap(pl4, metadata_subset,width=0.2,color = NULL,colnames_angle=45) +
    scale_fill_viridis_d(option="D")
#ggsave('figs/test_tree_iqtree_antibiotic.png', plot = pl5)







#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################



##############################
#markers ampicillin

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                #filter(CHROM == 'LMHECDEF_04343') %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "MFCAOJAD_04227_219" | snp_id == "MFCAOJAD_04226_396")

geno_rr_subset <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(MFCAOJAD_04227_219 = as.factor(MFCAOJAD_04227_219),
           MFCAOJAD_04226_396 = as.factor(MFCAOJAD_04226_396))


metadata_subset <- metadata %>% select( taxa, ampicillin) %>%
        mutate(ampicillin = as.factor(ampicillin))



geno_rr_subset2 <- geno_rr_subset %>% left_join(.,metadata_subset) %>% select(-taxa) %>% as.data.frame()

rownames(geno_rr_subset2) <- geno_rr_subset$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, geno_rr_subset2,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none")
#ggsave('figs/test_tree_iqtree_marker_ampicillin.png', plot = pl7)



#plot effects of markers
plot_markers <- geno_rr_subset2 %>% mutate(Resistance = ifelse(ampicillin == 'Resistant',1,
                                                ifelse(ampicillin == 'Intermediate',0.5,
                                                ifelse(ampicillin == 'Susceptible',0,NA)))) %>%
                                                        #mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(MFCAOJAD_04227_219,MFCAOJAD_04226_396))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw()

ggsave('figs/histogram_ampicillin_markers.png', plot = pl_effect, width = 5, height = 5)








#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


##############################
#markers trimethoprim.sulfamethoxazole

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                #filter(CHROM == 'LMHECDEF_04343') %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22" | snp_id == "DHJNCGMO_04404_672")

geno_rr_subset <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22` = as.factor(`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22`),
           DHJNCGMO_04404_672 = as.factor(DHJNCGMO_04404_672))


metadata_subset <- metadata %>% select( taxa, trimethoprim.sulfamethoxazole) %>%
        mutate(trimethoprim.sulfamethoxazole = as.factor(trimethoprim.sulfamethoxazole))



geno_rr_subset2 <- geno_rr_subset %>% left_join(.,metadata_subset) %>% select(-taxa) %>% as.data.frame()
table(geno_rr_subset2)
table(geno_rr_subset2$trimethoprim.sulfamethoxazole,geno_rr_subset2[,2])


rownames(geno_rr_subset2) <- geno_rr_subset$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, geno_rr_subset2,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none")
ggsave('figs/test_tree_iqtree_marker_trimethoprim.sulfamethoxazole.png', plot = pl7)



#plot effects of markers
plot_markers <- geno_rr_subset2 %>% mutate(Resistance = ifelse(trimethoprim.sulfamethoxazole == 'Resistant',1,
                                                ifelse(trimethoprim.sulfamethoxazole == 'Intermediate',0.5,
                                                ifelse(trimethoprim.sulfamethoxazole == 'Susceptible',0,NA)))) %>%
                                                        #mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(`Cluster_7017_+_+_GCA_003333865.1_ASM333386v1_+_+_DHJNCGMO_04403_+_+_DHJNCGMO_04404_+_+_DP_22`,DHJNCGMO_04404_672))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw()

#ggsave('figs/histogram_trimethoprim.sulfamethoxazole_markers.png', plot = pl_effect, width = 5, height = 5)geno_rr_subset2