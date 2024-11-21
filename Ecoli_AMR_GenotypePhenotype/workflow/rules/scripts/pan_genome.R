library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#df <- read.csv('pangenome/whole_pan_ecor_presence_absence.csv') %>% select(-X,-Group_annotation)
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]


df <- read.csv(input_file)

df2 <- df %>% select(-X)%>%
    pivot_longer(cols = c(-1)) %>% group_by(Locus, value) %>% summarize(count = n()) %>% filter(value == 'Present')


good_contigs <- df2 %>% filter(value == 'Present' & count == 72) %>% select(Locus) %>% unique()

#write.table(good_contigs,'pangenome/ecor_shared_contigs.txt', sep = '\t', row.names =F, col.names = F, quote=F)
#output list of good contigs to use in downstream analyses
write.table(good_contigs,output_file, sep = '\t', row.names =F, col.names = F, quote=F)
