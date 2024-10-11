#!/usr/bin/Rscript

library(SNPRelate)
library(ggplot2)
library(dplyr)
library(data.table)
library(gridExtra)
library(tidyr)

#input_vcf <- "vcf_files/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAF10_refedit_subsample.vcf"
#input_metadata <- "strains_metadata_phenotypes_full.txt"


###########################
input_vcf <- snakemake@input[['popgen_syn_vcf']]
input_metadata <- snakemake@input[['metadata_formatted']]

pca_fig_file <- snakemake@output[['pca_fig']]
pca_model_output_file <- snakemake@output[['pca_model_output']]




###########################
#create temporary gds file
snpgdsVCF2GDS(input_vcf, "test.gds", method="biallelic.only")

#read temp gds file
genofile <- snpgdsOpen("test.gds", readonly = F)

#read in metadata to explore which metadata features align with
metadata <- fread(input_metadata)




###########################
#run pca on variants
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

#extract eigencectors
tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2],  # the second eigenvector etc.
                    EV3 = pca$eigenvect[,3],
                    EV4 = pca$eigenvect[,4],
                    EV5 = pca$eigenvect[,5],
                    stringsAsFactors = FALSE)


#join PCA analysis results to metdata features of interest
#calculate new mlst (phylogroup) metdata column that only keeps an mlst label if there are enough observations (50)
df <- metadata %>% select(sample.id, Isolation.Country, Collection.Year,ref, Other.Typing) %>%
        right_join(.,tab, by = 'sample.id') %>%
                group_by(Other.Typing) %>% mutate(ntype = n()) %>%
                mutate(mlst = ifelse(ntype > 49, Other.Typing, NA) ) %>% ungroup(.)





###########################
#visualize PCA space with different metadata labels
pl_ref <- ggplot(df, aes(x=EV1, y=EV2, colour = ref)) + geom_point()
pl_mlst <-ggplot(df, aes(x=EV1, y=EV2, colour = mlst)) + geom_point()+ stat_ellipse()
pl_year <-ggplot(df, aes(x=EV1, y=EV2, colour = as.character(Collection.Year))) + geom_point() + labs(colour='Year')

#combine plots into one plot
pl2 <- grid.arrange(pl_ref, pl_year, pl_mlst, ncol=3, nrow =1)

#output figure
ggsave(pca_fig_file, pl2, width = 21)






###########################
#multinomial model on metadata features vs PCA
library(nnet)
library(pscl) #for pseudo R2, look at r2ML Maximum likelihood pseudo r-squared

#model of country ~ PCA
df$Isolation.Country<- as.factor(df$Isolation.Country)
mod1 <- multinom(Isolation.Country ~ EV1+EV2+EV3+EV4+EV5 , data = df,na.action=na.omit)
mod1r2 <- data.frame(t(pscl::pR2(mod1)))
mod1r2$model <- 'Isolation_country'

#model of collection year ~ PCA
df$Collection.Year <- as.factor(df$Collection.Year)
mod2 <- multinom(Collection.Year ~ EV1+EV2+EV3+EV4+EV5 , data = df,na.action=na.omit)
pscl::pR2(mod2)
mod2r2 <- data.frame(t(pscl::pR2(mod2)))
mod2r2$model <- 'Collection_year'

#model of mlst ~ PCA
df$mlst <- as.factor(df$mlst)
mod3 <- multinom(mlst ~ EV1+EV2+EV3+EV4+EV5 , data = df,na.action=na.omit)
pscl::pR2(mod3)
mod3r2 <- data.frame(t(pscl::pR2(mod3)))
mod3r2$model <- 'MLST'

#demonstrate that isolation by distance does emerge if we condition model on E. coli clade
mod4 <- multinom(Isolation.Country ~ EV1+EV2+EV3+EV4+EV5+mlst , data = df,na.action=na.omit)
pscl::pR2(mod4)
mod4r2 <- data.frame(t(pscl::pR2(mod4)))
mod4r2$model <- 'Isolation_country_conditioning_MLST'

#combine r2 estimates from all models
all_model_r2 <- rbind(mod1r2, mod2r2, mod3r2, mod4r2)

#output model fits to table
write.table(all_model_r2, pca_model_output_file, col.names = T, row.names = T, quote =F, sep = '\t' )