library(ggplot2)
library(dplyr)


df <- read.table("vcf_files/annotated_output_biallelic_SFS_counts.txt")


df <- df %>% mutate(annot = ifelse(grepl('HIGH', V5), 'LOF',
                                ifelse(grepl('synonymous', V5), 'synonymous',
                                    ifelse(grepl('missense', V5), 'missense',
                'NA'))))


sfs <- df %>% group_by(annot) %>%
        mutate(countannot = n()) %>%
            group_by(annot, V4) %>%
            summarize(AF = n()/countannot,
                      AC = n())




pl1 <- ggplot(sfs, aes(x=V4, y=AF, colour = annot)) + geom_point()

ggsave('SFS.pdf', pl1)



pl2 <- ggplot(sfs, aes(x=V4, y=AF, colour = annot)) + geom_smooth(method = "loess")
ggsave('SFS_loess.png', pl2)
