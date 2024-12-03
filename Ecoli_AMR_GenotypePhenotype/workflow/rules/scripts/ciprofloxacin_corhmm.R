install.packages("corHMM", version='2.8',repos = "http://cran.us.r-project.org")

library(ggtree)
library(dplyr)
library(data.table)
library(tidyr)
library(treeio)
library(corHMM)


#########
#snakemake file args
time_calibrated_tree_file <- snakemake@input[["gyrA_tree_time_calibrated"]]
marker_genotypes_file <- snakemake@input[["marker_genotypes"]]

output_model_file <- snakemake@output[["corhmm_ciprofloxaxin_markers"]]


#########
#read in data, time calibrated species tree and marker state for top SNP markers
time_calibrated_tree <- read.newick (time_calibrated_tree_file)
geno <- fread(marker_genotypes_file)





############################################################
############################################################
#markers

#funs to clean bcftools column names and convert missing genotyprs to reference
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }

#converts missing genotypes to reference
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
                 filter(snp_id == "LMHECDEF_04343_259" | snp_id == "LMHECDEF_04343_248"|  snp_id == "LMHECDEF_01124_239") #selects focal markers

#convert to long format, code genotypes as factor and rename markers to gene name + position
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

#convert missing genotypes to reference
geno_cip_markers_long <- convert_missing_ref_GT(geno_cip_markers_long)


#####################################################################################
#####################################################################################


#prepare for input to corHMM< for gyrA259 where there is 3 alternate alleles, recode as all alternates are '1'
hmm_data <- geno_cip_markers_long%>%
    mutate(gyrA_259 = ifelse(as.numeric(as.character(gyrA_259)) > 0,1,0)) %>%
    mutate(gyrA_259 = as.factor(gyrA_259),
           gyrA_248 = as.factor(gyrA_248),
           parC_239 = as.factor(parC_239)) %>% data.frame(.)

#run corHMM with 1 rate category
mk_mod <- corHMM(phy = time_calibrated_tree, data = hmm_data, rate.cat = 1)
mk_mod

#save model output
saveRDS(mk_mod, output_model_file)
