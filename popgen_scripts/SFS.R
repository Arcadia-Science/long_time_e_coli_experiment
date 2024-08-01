library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)

#input counts scarped from VCF (assuming missing GT is ref)
df <- read.table("vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts.txt")

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
                      AC = n())



#plot SFS
pl1 <- ggplot(sfs, aes(x=V4, y=AF, colour = annot)) + geom_point() + geom_line()  + xlim(0,20)

ggsave('SFS_goodcontigs.pdf', pl1)

############################
##look at coverage by dividing allele count by total depth

dpth <- df %>% mutate(depth = V5/V4)

#all sites
pl_depth <- ggplot(dpth, aes(x=depth, colour = annot)) + geom_density()
pl_count_Depth <- ggplot(dpth, aes(x=V4, y=V5, colour = annot)) + geom_point(alpha = 0.1 )
pl2 <- grid.arrange(pl_depth,pl_count_Depth, ncol = 2)

ggsave('Seq_Depth_figs.png', pl2, width = 22, height = 9)

#focus on rare
pl_depth_rare <- ggplot(dpth, aes(x=depth, colour = annot)) + geom_density()
pl_count_Depth_rare <- ggplot(dpth, aes(x=V4, y=V5, colour = annot)) + geom_point(alpha = 0.1 ) + xlim(0,5)
pl2_rare <- grid.arrange(pl_depth_rare,pl_count_Depth_rare, ncol = 2)

ggsave('Seq_Depth_figs_rare.png', pl2_rare, width = 22, height = 9)
###############################

##look closer at quadrupletons

dpth2 <- dpth %>% filter(annot == 'synonymous' & V4 %in% c(1:4)) %>%
            group_by(V1, V4) %>%
            summarize(count = n(),
                      contig_Depth = mean(depth, na.rm=T)) %>%
                        pivot_wider(names_from = V4, values_from = count)

#calculate ratio of site types to quadrupletons per contig

dpth2 <- dpth2 %>% mutate(rat1 = `4`/`1`,
                         rat2 = `4`/`2`,
                         rat3 = `4`/`3`)

pl_rat1 <- ggplot(dpth2, aes(x=rat1)) + geom_density()
pl_quad_contig_depth <- ggplot(dpth2, aes(x=contig_Depth, y=rat1)) + geom_point( )

pl3 <- grid.arrange(pl_rat1,pl_quad_contig_depth, ncol = 2)

ggsave('quadrupletons_ratio_per_contig.png', pl3, width = 10, height = 7)
