library(data.table)
library(dplyr)
library(tidyr)
library(rrBLUP)


#head(geno2[,1:7])


#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
#convert genotype values for rrBLUP to conform to -1 = 0/0
convert_gt_value <- function(x){ gsub("0", "-1", x)}


#input data
pheno <- fread('strains_metadata_phenotypes_full.txt')

#geno <- fread('vcf_files/annotated_output_biallelic_synonymous_MAF7_refedit_goodcontigs_subsample_geno.txt')
geno <- fread('geno_pred/annotated_output_biallelic_refedit_goodcontigs_ldpruned_remoutliers_edited_subsample.txt')




######################################################
#prep pheno

pheno2 <- pheno %>% select(sample.id, ciprofloxacin)



######################################################
#prep geno


#clean column names from bcftools output and create relevant columns
geno_rr <- geno %>%
                mutate(snp_id = paste(`#CHROM`, POS, sep = '_')) %>%
                relocate(snp_id, .before = `#CHROM`) %>%
                select(-`#CHROM`, -POS, -ALT, -AC, -DP)



#tranpose to proper orientation
geno_rr_2 <- geno_rr%>%
  pivot_longer(cols=c(-snp_id),names_to="sample.id")%>%
  pivot_wider(names_from=c(snp_id))

#head(geno_rr_2[,1:7])

######################################################
#head(pheno_geno[,1:7])
pheno_geno <- inner_join(pheno2, geno_rr_2) %>% filter(!is.na(ciprofloxacin))

input_pheno <- pheno_geno %>% select(ciprofloxacin) %>% mutate(phen = ifelse(ciprofloxacin == 'Resistant',1,0)) %>% select(phen) #%>% as.matrix()
input_pheno2 <- as.matrix(input_pheno)


input_geno <- pheno_geno %>% select(-sample.id, -ciprofloxacin) %>% mutate_if(is.character, as.numeric) #%>% as.matrix(.)
input_geno2 <- as.matrix(input_geno)

model_rrblup <- mixed.solve(input_pheno2, Z = input_geno2, K=NULL, SE = FALSE, return.Hinv=FALSE)

snp_effects <- model_rrblup$u

######################################################
df <- as.data.frame(snp_effects) %>% arrange(desc(snp_effects))

df <- as.data.frame(snp_effects) %>% mutate(chrom = gsub('_([^_]*)$','',rownames(.)),
                                            pos = gsub('.*\\_','',rownames(.)))%>%
                                              arrange(desc(snp_effects)) %>% left_join(.,gene_names)


df2 <- df %>% group_by(chrom) %>% summarize(count_snp = n(), effect = mean(snp_effects)) %>% arrange(desc(effect))

gene_names <- fread('pangenome/gene_names.txt', fill=TRUE, header = F) %>% rename(chrom = V1, gene = V2)

df