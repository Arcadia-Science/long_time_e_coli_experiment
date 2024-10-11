library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggnewscale)
library(ape)
library(castor)
library(geiger)

#tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.treefile")
tree_test <- read.newick("tree/alignment_remoutliers_syn_MAC10/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.min4.phy.timetree.nwk")

tree_poly_resolved <- multi2di(tree_test)

metadata <- fread("strains_metadata_phenotypes_full.txt")
geno <- fread('geno_pred/annotated_output_MAF250_gwas_bslmm_probit_top_hits_gts.txt')

metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )


metadata_subset <- metadata %>% select(taxa, ciprofloxacin, ampicillin, trimethoprim.sulfamethoxazole) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin),
                ampicillin = as.factor(ampicillin),
                trimethoprim.sulfamethoxazole = as.factor(trimethoprim.sulfamethoxazole))


dtt_data <- metadata_subset %>% select(ciprofloxacin)%>%
        mutate(testpheno = ifelse(ciprofloxacin == 'Susceptible' | is.na(ciprofloxacin),0,1)) %>%
                select(-ciprofloxacin) %>%
                as.matrix(.)
rownames(dtt_data) <- metadata_subset$taxa

dtt_antibiotics = dtt(phy = tree_poly_resolved, data=dtt_data)


#######ltt and dtt analyses
#lineage through time plot
ltt.plot(tree_test)

pdf(file="figs/ltt_phylo_mac10.pdf")
ltt.plot(tree_test)
dev.off()

#disparity through time plots



trees<-pbtree(n=100,scale=1,nsim=1)
