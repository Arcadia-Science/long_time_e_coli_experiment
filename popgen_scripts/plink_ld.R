library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)


df <- fread('plink/annotated_output_biallelic_synonymous_MAF7_refedit_subsample.ld')
df$distance <- abs(df$BP_A - df$BP_B )



pl <- ggplot(df, aes(x=distance, y=R2)) + geom_point(alpha=0.1) + geom_smooth(method = "loess", se=FALSE)
ggsave('ld_decay.png', pl)



pl2 <- ggplot(df, aes(x=distance, y=R2)) + geom_point(alpha=0.1) + geom_smooth(method = "loess", se=FALSE) + xlim(0,1000)
ggsave('ld_decay_1kb.png', pl2)