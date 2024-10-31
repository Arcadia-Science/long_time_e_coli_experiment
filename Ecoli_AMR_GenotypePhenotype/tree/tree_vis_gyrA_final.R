library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggnewscale)
library(corHMM)

gyrA_tree<- read.newick ("tree/alignment_gyra/annotated_output_gyrA_refedit.min4.phy.treefile")

metadata <- fread("strains_metadata_phenotypes_full.txt")
geno <- fread('geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt')

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )




############################################################
############################################################
#markers

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#prepare metadata and marker genotypes for visualization
geno_cip_markers <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP) %>%
                 filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248"|  snp_id == "LMHECDEF_01124_239")

geno_cip_markers_long <- geno_cip_markers%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(LMHECDEF_04343_259 = as.factor(LMHECDEF_04343_259),
           LMHECDEF_04343_248 = as.factor(LMHECDEF_04343_248),
           LMHECDEF_01124_239 = as.factor(LMHECDEF_01124_239)) %>%
                rename(
                        gyrA_259 = LMHECDEF_04343_259,
                        gyrA_248 = LMHECDEF_04343_248,
                        parC_239 = LMHECDEF_01124_239)


pheno_cip <- metadata %>% select( taxa, ciprofloxacin) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin))



cip_metdata_full <- geno_cip_markers_long %>% left_join(.,pheno_cip) %>% select(-taxa) %>% as.data.frame()
rownames(cip_metdata_full) <- geno_cip_markers_long$taxa


######################
#plot tree and metadata



pl_base <- ggtree(gyrA_tree,layout='circular')


gyra_tree_plot <- gheatmap(pl_base, cip_metdata_full,width=0.2,color = NULL,colnames_angle=45) +
 scale_fill_manual(breaks=c("0", "1","2","3", "Intermediate", "Resistant", "Susceptible"),
        values=c("#73B5E3", "#F7B846","#3B9886","#97CD78" , "#B5BEA4", "#F28360", "#5088C5"),  na.value = NA)


ggsave('final_figs/Fig5A_gyrA_tree_500.png', plot = gyra_tree_plot)
ggsave('final_figs/Fig5A_gyrA_tree_500.svg', plot = gyra_tree_plot)








################################################################
################################################################

#plot effects of markers gyrA



#plot effects of top 3 markers
plot_markers <- cip_metdata_full %>% mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                        mutate(gyrA_259_simp = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(gyrA_248,gyrA_259_simp,parC_239))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        xlab('Genotype')+
        scale_x_discrete(labels=c(
                                "000" = "Ancetral\nstate",
                                "001" = "parC239",
                                "010" = "gyrA259",
                                "100" = "gyrA248",
                                "101" = "gyrA248\nparC239",
                                "110" = "gyrA248\ngyrA259",
                                "111" = "Full stack"
                                ))

ggsave('final_figs/Fig5B_boxplot_marker_effects_500.png', plot = pl_effect, width = 7, height = 5)
ggsave('final_figs/Fig5B_boxplot_marker_effects_500.svg', plot = pl_effect, width = 7, height = 5)




#############
epi <- cip_metdata_full %>% mutate(pheno = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA))) %>%
                                                mutate(gyrA_259_simp = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
                                                mutate( gyrA_248 = as.numeric(as.character(gyrA_248)),
                                                        gyrA_259_simp = as.numeric(as.character(gyrA_259_simp)),
                                                        parC_239 = as.numeric(as.character(parC_239)))

mylogit <- glm(pheno ~ gyrA_248 + gyrA_259_simp + parC_239 + gyrA_248*gyrA_259_simp + parC_239*gyrA_259_simp + parC_239*gyrA_248,
                data = epi, family = 'binomial')

summary(mylogit)





################################################################
################################################################
#corHMM state transition estimation for ciprofloxacin top 3 markers

#full_tree <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")
full_tree <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.timetree.nwk")


hmm_data <- geno_cip_markers_long%>%
    mutate(gyrA_259 = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
    mutate(gyrA_259 = as.factor(gyrA_259),
           gyrA_248 = as.factor(gyrA_248),
           parC_239 = as.factor(parC_239)) %>% data.frame(.)


mk_mod <- corHMM(phy = full_tree, data = hmm_data, rate.cat = 1)
mk_mod


saveRDS(mk_mod, file = "tree/hmm_ciprofloxacin_fulltree_final.rda")

######graph plotting
library(igraph)

rates <- as.matrix(mk_mod$solution)
colnames(rates) <- c("WT","parE249", "gyrA248", "gyrA248_parE249", "gyrA259", "gyrA248_gyrA259", "full_stack")
colnames(rates) <- c("WT","parE249", "gyrA248", "gyrA248_parE249", "gyrA259", "gyrA248_gyrA259", "full_stack")

#colnames(rates) <- c("0_0_0","0_0_1", "0_1_0", "0_1_1", "1_0_0", "1_1_0", "1_1_1")
#rownames(rates) <- c("0_0_0","0_0_1", "0_1_0", "0_1_1", "1_0_0", "1_1_0", "1_1_1")

rates[2,1] <- 0
rates[3,1] <- 0
rates[5,1] <- 0

rates[4,2] <- 0
rates[4,3] <- 0
rates[6,3] <- 0
rates[6,5] <- 0

rates[7,4] <- 0
rates[7,6] <- 0

rates[is.na(rates)] <- 0


testgraph <- graph_from_adjacency_matrix(
  rates,
  mode = c("directed"),
  weighted = TRUE,
  diag = F,
  add.colnames = NULL,
  add.rownames = NA
)

mylay <- layout_on_grid(testgraph)
customlay = c(
0,1,
1,0,
1,1,
2,0.5,
1,2,
2,1.5,
3,1)

dim(customlay) <- c(2, 7)
customlay <- t(customlay)







E(testgraph)$width <- E(testgraph)$weight + min(E(testgraph)$weight) + 1
pl_graph <- plot(testgraph)
#ggsave('figs/TESTGRAPH.png', plot = pl_graph)

pdf(file="figs/TESTGRAPH.pdf")
plot(testgraph, edge.width=log(E(testgraph)$weight), edge.curved = 0.3, layout=customlay)
dev.off()
