#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)


#snakemake file arguments
input_sfs_file <- snakemake@input[['sfs_no_outlier_input']]
input_presence_absence_sfs_file <- snakemake@input[['presence_absence_input']]
input_presence_absence_sfs_contigs <- snakemake@input[['presence_absence_loci']]


output_sfs_file <- snakemake@input[['sfs_fig']]


#########################
#bi allelic SNP SFS

#df <- read.table('vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts_remoutliers.txt')


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
pl_sfs <- ggplot(sfs, aes(x=V4, y=AF, colour = annot)) +
    geom_point() + geom_line() + xlim(0,20) + theme_bw() +
    xlab('Allele Count') + ylab('Proportion of SNPs')




#########################
#presence_absence SFS

#sfs <- fread('presence_absence_data/presence_absence_sfs.txt')%>% select(-contains("V"))
#chroms <- fread('presence_absence_data/presence_absence_all.bim') %>% select(V1)


#read in table of presence absence calls, and name of contigs (not strictly needed for plot but nice to have)
sfs <- fread(input_presence_absence_sfs_file)%>% select(-contains("V"))
chroms <- fread(input_presence_absence_sfs_contigs) %>% select(V1)

#merge and calculate counts across all individuals
df <- cbind(chroms, sfs)
df$present <-  rowSums(df[,2:7044] == 2,na.rm=T)
df$absent <-  rowSums(df[,2:7044] == 1,na.rm=T)

#subset to new clean dataframe and check that the totals all add up to 7043
df2 <- df %>% select(V1, present, absent)
df2$total <- df2$present + df2$absent
table(df2$total)

#plot presence absence
pl_pres_abs <- ggplot(df2, aes(x=present)) + geom_histogram() + theme_bw()




#########################
#save plot

pl_combined <- grid.arrange(pl_sfs,  pl_pres_abs, ncol=2, nrow =1)

ggsave(output_sfs_file,pl_combined, height = 4, width = 10.5)
#ggsave('figs/SFS_combined_remoutliers.png',pl_combined, height = 4, width = 10.5)
