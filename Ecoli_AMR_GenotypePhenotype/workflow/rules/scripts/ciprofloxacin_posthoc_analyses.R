library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggnewscale)




#gyrA_tree<- read.newick ("../tree/alignment/gyrA_remoutliers_edited.min4.phy.treefile")
#metadata <- fread("../strains_metadata_phenotypes_full.txt")
#geno <- fread('../geno_pred/gwas_top10_marker_genotypes.txt')

#########
#snakemake file args
gyrA_tree_file <- snakemake@input[["gyrA_tree"]]
metadata_formatted_file <- snakemake@input[["metadata_formatted"]]
marker_genotypes_file <- snakemake@input[["marker_genotypes"]]

gyrA_tree_fig <- snakemake@output[["gyrA_tree_fig"]]
cipro_marker_boxplot <- snakemake@output[["cipro_marker_boxplot"]]
cipro_marker_logistic_regression <- snakemake@output[["cipro_marker_logistic_regression"]]



#########
#read in data
gyrA_tree<- read.newick (gyrA_tree_file)
metadata <- fread(metadata_formatted_file)
geno <- fread(marker_genotypes_file)





###############
metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )




############################################################
############################################################
#markers

#funs to cleat bcftools column names and convert missing genotyprs to reference
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }

convert_missing_ref_GT  <- function(df) {
  # Apply gsub over each cell in the data frame
  df[] <- lapply(df, function(column) {
    # Check if the column is character or factor, then replace "." with "0"
    if (is.character(column) || is.factor(column)) {
      gsub("\\.", "0", as.character(column))
    } else {
      column  # Return the column unchanged if it's not character or factor
    }
  })
  return(df)
}



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
                        parC_239 = LMHECDEF_01124_239) %>% convert_missing_ref_GT(.)


geno_cip_markers_long <- convert_missing_ref_GT(geno_cip_markers_long)


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


ggsave(gyrA_tree_fig, plot = gyra_tree_plot)








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

ggsave(cipro_marker_boxplot, plot = pl_effect, width = 7, height = 5)




#############
epi <- cip_metdata_full %>% mutate(pheno = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA))) %>%
                                                mutate(gyrA_259_simp = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
                                                mutate( gyrA_248 = as.numeric(as.character(gyrA_248)),
                                                        gyrA_259_simp = as.numeric(as.character(gyrA_259_simp)),
                                                        parC_239 = as.numeric(as.character(parC_239)))

mylogit <- glm(pheno ~ gyrA_248 + gyrA_259_simp + parC_239 + gyrA_248*gyrA_259_simp + parC_239*gyrA_259_simp + parC_239*gyrA_248,
                data = epi, family = 'binomial')

modsum <- summary(mylogit)$coefficients


write.table(modsum, cipro_marker_logistic_regression)
