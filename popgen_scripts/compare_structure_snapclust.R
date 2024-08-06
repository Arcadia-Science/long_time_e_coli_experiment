library(data.table)
library(dplyr)
library(ggplot2)

df1 <- fread('tables/snapclust_synonymous_goodcontigs_k10_results.txt')

df2 <- fread('tables/structure_synonymous_goodcontigs_k9_k10_results.txt')


df<- full_join(df1, df2, by = c('K','sampID','ID','Season','isolation_country','collection_year')) %>% filter(!is.na(Cluster.y) & !is.na(Cluster.x))


table(df$Cluster.y, df$Cluster.x)
