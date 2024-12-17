library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

df <- read.csv(input_file)
number_ECOR_strains <- 72

df2 <- df %>% select(-X) %>%
    pivot_longer(cols = c(-1)) %>% group_by(Locus, value) %>% summarize(count = n()) %>% filter(value == 'Present')

#check that number of obersvations matches number of strains in dataset (72) for Presence values
good_contigs <- df2 %>% filter(value == 'Present' & count == number_ECOR_strains) %>% select(Locus) %>% unique()

#output list of good contigs to use in downstream analyses
write.table(good_contigs, output_file, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
