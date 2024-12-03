library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggnewscale)





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
gyrA_tree<- read.newick (gyrA_tree_file) #gene tree
metadata <- fread(metadata_formatted_file) #metadata
geno <- fread(marker_genotypes_file) #marker genotypes





###############
metadata <- metadata %>%
        rename(taxa = sample.id) %>%
        relocate(taxa, .before = X )




############################################################
############################################################
#markers

#funs to clean bcftools column names and convert missing genotyprs to reference
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }

#convert missing genotypes to reference allele
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
                 filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248"|  snp_id == "LMHECDEF_01124_239") #select focal markers

#convert to long and recode as factor
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

#convert missing genotypes to reference allele
geno_cip_markers_long <- convert_missing_ref_GT(geno_cip_markers_long)

#select sample names (to add to tree tips) and ciprofloxacin phenotype, recode as factor
pheno_cip <- metadata %>% select( taxa, ciprofloxacin) %>%
        mutate(ciprofloxacin = as.factor(ciprofloxacin))


#add ciprofloxacin phenotype to focal markers, rename rows to sample names to plot on tree
cip_metdata_full <- geno_cip_markers_long %>% left_join(.,pheno_cip) %>% select(-taxa) %>% as.data.frame()
rownames(cip_metdata_full) <- geno_cip_markers_long$taxa


######################
#plot tree and metadata


#base circular gene tree
pl_base <- ggtree(gyrA_tree,layout='circular')

#plot tree with ciprofloxacin phenotype and 3 focal markers (allele state) added on
gyra_tree_plot <- gheatmap(pl_base, cip_metdata_full,width=0.2,color = NULL,colnames_angle=45) +
 scale_fill_manual(breaks=c("0", "1","2","3", "Intermediate", "Resistant", "Susceptible"),
        values=c("#73B5E3", "#F7B846","#3B9886","#97CD78" , "#B5BEA4", "#F28360", "#5088C5"),  na.value = NA)


ggsave(gyrA_tree_fig, plot = gyra_tree_plot)








################################################################
################################################################

#plot effects of markers gyrA



#plot boxplot of effects of top 3 markers
#code resistance as 0=susceptible, 0.5 = intermediate, and 1=resistant
#code all alternate alleles of gyrA259 to be equivalent
plot_markers <- cip_metdata_full %>% mutate(Resistance = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Intermediate',0.5,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA)))) %>%
                                                        mutate(gyrA_259_simp = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
                                                                mutate(genotype = paste0(gyrA_248,gyrA_259_simp,parC_239))

#boxplot of ciprofloxacin resistance phenotype vs allele state (ref/alt)
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
#prepare input for logistic regression model to test for interactions between markers and ciprofloxacin resistance
#code resistance as 0=susceptible, and 1=resistant, ignore rare intermediate (NA)
#code all alternate alleles of gyrA259 to be equivalent
epi <- cip_metdata_full %>% mutate(pheno = ifelse(ciprofloxacin == 'Resistant',1,
                                                ifelse(ciprofloxacin == 'Susceptible',0,NA))) %>%
                                                mutate(gyrA_259_simp = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
                                                mutate( gyrA_248 = as.numeric(as.character(gyrA_248)),
                                                        gyrA_259_simp = as.numeric(as.character(gyrA_259_simp)),
                                                        parC_239 = as.numeric(as.character(parC_239)))

#run logistic regression of allele state of arkers vs ciprofloxacin resistance include all pairwise interaction terms
mylogit <- glm(pheno ~ gyrA_248 + gyrA_259_simp + parC_239 + gyrA_248*gyrA_259_simp + parC_239*gyrA_259_simp + parC_239*gyrA_248,
                data = epi, family = 'binomial')

#output summary of model to table
modsum <- summary(mylogit)$coefficients
write.table(modsum, cipro_marker_logistic_regression)
