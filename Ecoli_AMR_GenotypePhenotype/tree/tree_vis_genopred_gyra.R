library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggnewscale)
library(corHMM)

tree_test<- read.newick ("tree/alignment_gyra/annotated_output_gyrA_refedit.min4.phy.treefile")

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
#ggsave('figs/test_tree_iqtree_gyra_gamma.png', plot = pl)


pl4 <- pl + geom_tippoint(aes(colour = mlst),size = 0.001) +  theme(legend.position = "none")
#ggsave('figs/test_tree_iqtree_gyra_mlst.png', plot = pl4)
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


############################################################
############################################################
#markers

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                filter(CHROM == 'LMHECDEF_04343' | CHROM ==  "LMHECDEF_01124") %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248"|  snp_id == "LMHECDEF_01124_239")

geno_rr_subset <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(LMHECDEF_04343_259 = as.factor(LMHECDEF_04343_259),
           LMHECDEF_04343_248 = as.factor(LMHECDEF_04343_248),
           LMHECDEF_01124_239 = as.factor(LMHECDEF_01124_239))


metadata_subset <- metadata %>% select( taxa, ciprofloxacin) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin))



geno_rr_subset2 <- geno_rr_subset %>% left_join(.,metadata_subset) %>% select(-taxa) %>% as.data.frame()

rownames(geno_rr_subset2) <- geno_rr_subset$taxa
pl4 <- ggtree(tree_test,layout='circular')

pl7 <- gheatmap(pl4, geno_rr_subset2,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none")
ggsave('figs/test_tree_iqtree_gyra_3markers.png', plot = pl7)


#### custom colours
pl7 <- gheatmap(pl4, geno_rr_subset2,width=0.2,color = NULL,colnames_angle=45) + theme(legend.position = "none") +
 scale_fill_manual(breaks=c("0", "1","2","3", "Intermediate", "Resistant", "Susceptible"),
        values=c("steelblue", "darkorange3","red2","red3" , "red", "goldenrod3", "skyblue"),  na.value = NA)


ggsave('figs/test_tree_iqtree_gyra_3markers.png', plot = pl7)





################################################################
################################################################

#plot effects of markers gyrA
plot_markers <- geno_rr_subset2 %>% mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                        mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(LMHECDEF_04343_248,LMHECDEF_04343_259))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw()

ggsave('figs/histogram_ciprofloxacin_markers.png', plot = pl_effect, width = 5, height = 5)





#plot effects of top 3 markers
plot_markers <- geno_rr_subset2 %>% mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                        mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(LMHECDEF_04343_248,LMHECDEF_04343_259,LMHECDEF_01124_239))

pl_effect <- ggplot(plot_markers, aes(x=genotype, y = Resistance)) + geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.05)) + theme_bw()

ggsave('figs/histogram_ciprofloxacin_markers_top3.png', plot = pl_effect, width = 5, height = 5)


#############
epi <- geno_rr_subset2 %>% mutate(pheno = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA))) %>%
                                                mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
                                                mutate( gyrA_248 = as.numeric(as.character(LMHECDEF_04343_248)),
                                                        gyrA_259 = as.numeric(as.character(LMHECDEF_04343_259)),
                                                        parC_239 = as.numeric(as.character(LMHECDEF_01124_239)))

mylogit <- glm(pheno ~ gyrA_248 + gyrA_259 + parC_239 + gyrA_248*gyrA_259 + parC_239*gyrA_259 + parC_239*gyrA_248,
                data = epi, family = 'binomial')

summary(mylogit)

fit_land <- epi %>% mutate(hamming_dist = gyrA_248+gyrA_259+parC_239) %>%
                        mutate(geno = paste0(gyrA_248,gyrA_259,parC_239)) %>%
                        mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                group_by(hamming_dist, geno) %>%
                                                summarize(avg_resistance = mean(Resistance, na.rm=T),
                                                        count_geno = n())


pl_fl <- ggplot(fit_land, aes(y=avg_resistance, x = hamming_dist, colour = geno)) + geom_point(aes(size = count_geno))  + theme_bw() +
        geom_segment(aes(x = 0,y = 0.00897,xend = 1,yend = 0.189),arrow=arrow(length=unit(.1,'cm')), lwd = 0.5, color = 'grey') +
        geom_segment(aes(x = 1,y = 0.189,xend = 2,yend = 0.955 ),arrow=arrow(length=unit(.1,'cm')), lwd = 0.5, color = 'grey') +
         geom_segment(aes(x = 2,y = 0.955,xend = 3,yend = 0.986 ),arrow=arrow(length=unit(.1,'cm')), lwd = 0.5, color = 'grey')

ggsave('figs/fit_land_ciprofloxacin_markers_top3.png', plot = pl_fl, width = 5, height = 5)


head(fit_land[fit_land$geno== '010',])

####################
#corHMM state transition estimation for ciprofloxacin top 3 markers

#full_tree <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")
full_tree <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.timetree.nwk")

geno <- fread('geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt')

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }


#marker genotypes in traget gene
geno_rr <- geno %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                filter(CHROM == 'LMHECDEF_04343' | CHROM ==  "LMHECDEF_01124") %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP)

geno_rr_targets <- geno_rr %>% filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248"|  snp_id == "LMHECDEF_01124_239")

hmm_data <- geno_rr_targets%>%
  group_by(snp_id) %>%
    pivot_longer(cols=c(-snp_id),names_to="taxa")%>%
    pivot_wider(names_from=c(snp_id)) %>%
    mutate(LMHECDEF_04343_259 = ifelse(as.numeric(as.character(LMHECDEF_04343_259)) > 0,1,0)) %>%
    mutate(LMHECDEF_04343_259 = as.factor(LMHECDEF_04343_259),
           LMHECDEF_04343_248 = as.factor(LMHECDEF_04343_248),
           LMHECDEF_01124_239 = as.factor(LMHECDEF_01124_239)) %>% data.frame(.)




mk_mod <- corHMM(phy = full_tree, data = hmm_data, rate.cat = 1)
mk_mod

pdf(file="figs/hmm_simplestates_ciprofloxacin_fulltree.pdf")
plotMKmodel(mk_mod)
dev.off()

saveRDS(mk_mod, file = "tree/hmm_simplestates_ciprofloxacin_fulltree.rda")

########
#forward rates only

custom_statemat <- getStateMat4Dat(hmm_data)
custom_statemat_forwardonly <- dropStateMatPars(custom_statemat$rate.mat, c(1,2,3,5,7,8,11,13,16))

#getFullMat(custom_statemat)

mk_mod_forward <- corHMM(phy = full_tree, data = hmm_data, rate.cat = 1, rate.mat = custom_statemat_forwardonly)
mk_mod_forward


custom_statemat_double <- getStateMat4Dat(hmm_data)
custom_statemat_double$rate.mat[6,7] <- 0
mk_mod_double <- corHMM(phy = full_tree, data = hmm_data, rate.cat = 1, rate.mat = custom_statemat_double)
#check_save <- readRDS("tree/hmm_allstates_ciprofloxacin_gyratree.rda")





###########################
#castor makrkov model rates
#prepare tree tip data

hmm_data2 <- hmm_data %>% mutate(geno = paste0(LMHECDEF_04343_248, LMHECDEF_04343_259, LMHECDEF_01124_239)) %>%
        mutate(geno_tips = as.numeric(as.factor(geno))) %>% select(taxa, geno_tips)

#states for ease of interpretability
states <- hmm_data %>% mutate(geno = paste0(LMHECDEF_04343_248, LMHECDEF_04343_259, LMHECDEF_01124_239)) %>%
        mutate(geno_tips = as.numeric(as.factor(geno))) %>% select(geno, geno_tips) %>% unique(.)

#order tree tips same as the phylogeny being used
tree_tips_order <- data.frame(full_tree$tip.label) %>% rename(taxa = full_tree.tip.label) %>%
                left_join(.,hmm_data2, by = 'taxa')

tree_tips_input <- c(tree_tips_order$geno_tips)







#one rate for all transitions
results_er = fit_mk(full_tree, 7, tree_tips_input, rate_model="ER",
                 Ntrials=5, optim_max_iterations = 5000, Nthreads = 20)

#any transition
results_ard = fit_mk(full_tree, 7, tree_tips_input, rate_model="ARD",
                 Ntrials=5, optim_max_iterations = 5000, Nthreads = 20)

#any transition but fix forward/reverse to be symmetrical
results_sym = fit_mk(full_tree, 7, tree_tips_input, rate_model="SYM", Ntrials=5,
                 optim_max_iterations = 5000, Nthreads = 20)


#only allow singleton jumps in hamming distance
single_only_mat <- results_sym$transition_matrix
single_only_mat[,1] <- c(0,1,2,3,0,0,0)
single_only_mat[,2] <- c(4,0,0,0,5,0,0)
single_only_mat[,3] <- c(6,0,0,0,0,7,0)
single_only_mat[,4] <- c(8,0,0,0,9,10,0)
single_only_mat[,5] <- c(0,11,0,12,0,0,13)
single_only_mat[,6] <- c(0,0,14,15,0,0,16)
single_only_mat[,7] <- c(0,0,0,0,17,18,0)

results_single_only = fit_mk(full_tree, 7, tree_tips_input, rate_model=single_only_mat,
                 Ntrials=5, optim_max_iterations = 5000, Nthreads = 20)


#only allow singleton jumps and only allow forward mutation 0>1
single_forward_only_mat <- results_sym$transition_matrix
single_forward_only_mat[,1] <- c(0,1,2,3,0,0,0)
single_forward_only_mat[,2] <- c(0,0,0,0,4,0,0)
single_forward_only_mat[,3] <- c(0,0,0,0,0,5,0)
single_forward_only_mat[,4] <- c(0,0,0,0,6,7,0)
single_forward_only_mat[,5] <- c(0,0,0,0,0,0,8)
single_forward_only_mat[,6] <- c(0,0,0,0,0,0,9)
single_forward_only_mat[,7] <- c(0,0,0,0,0,0,0)

results_single_forward_only = fit_mk(full_tree, 7, tree_tips_input, rate_model=single_forward_only_mat,
                 Ntrials=5, optim_max_iterations = 5000,, Nthreads = 20)



#allow singleton or doubleton jumps
double_single_mat <- results_sym$transition_matrix
double_single_mat[,1] <- c(0,1,2,3,4,5,0)
double_single_mat[,2] <- c(6,0,7,8,9,0,10)
double_single_mat[,3] <- c(11,12,0,13,0,14,15)
double_single_mat[,4] <- c(16,17,18,0,19,20,21)
double_single_mat[,5] <- c(22,23,0,24,0,25,26)
double_single_mat[,6] <- c(27,0,28,29,30,0,31)
double_single_mat[,7] <- c(0,32,33,34,35,36,0)

results_double_single = fit_mk(full_tree, 7, tree_tips_input, rate_model=double_single_mat,
                 Ntrials=5, optim_max_iterations = 5000, , Nthreads = 20)




model_comp <- data.frame(aic_mod = c(
        results_er$AIC,
        results_ard$AIC,
        results_sym$AIC,
        results_single_only$AIC,
        results_single_forward_only$AIC,
        results_double_single$AIC),

        n_iterations = c(
        results_er$Niterations,
        results_ard$Niterations,
        results_sym$Niterations,
        results_single_only$Niterations,
        results_single_forward_only$Niterations,
        results_double_single$Niterations),

        loglikelihood = c(
        results_er$loglikelihood,
        results_ard$loglikelihood,
        results_sym$loglikelihood,
        results_single_only$loglikelihood,
        results_single_forward_only$loglikelihood,
        results_double_single$loglikelihood),

                        mod_name = c(
                        'one_rate',
                        'all_free_rate',
                        'all_symmetrical',
                        'single_all',
                        'single_forward_only',
                        'results_double_single' )
)


#write.table(model_comp, 'tables/ciprofloxacin_top3_markers_markov_transitions.txt', col.names = T, quote = F, row.names = F, sep = '\t')






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
