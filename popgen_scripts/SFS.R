#!/usr/bin/Rscript


library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)

#input counts scarped from VCF (assuming missing GT is ref)
#df <- read.table("vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts_remoutliers.txt")

input_sfs_file <- snakemake@input[['sfs_no_outlier_input']]
output_sfs_file <- snakemake@input[['sfs_no_outlier_input']]

df <- read.table(input_sfs_file)


#simplify SNPeff annotation
df <- df %>% mutate(annot = ifelse(grepl('HIGH', V6), 'LOF',
                                ifelse(grepl('synonymous', V6), 'synonymous',
                                    ifelse(grepl('missense', V6), 'missense',
                'NA'))))

#calculate numbers for SFS bins and calculate frequency for each type of site
sfs <- df %>% group_by(annot) %>%
        mutate(countannot = n()) %>%
            group_by(annot, V4) %>%
            summarize(AF = n()/countannot,
                      AC = n()) %>% filter(V4 > 0)



#plot SFS
pl1 <- ggplot(sfs, aes(x=V4, y=AF, colour = annot)) +
    geom_point() + geom_line() + xlim(0,20) + theme_bw() +
    xlab('Allele Count') + ylab('Proportion of SNPs')

ggsave(output_sfs_file, pl1, height = 5)
