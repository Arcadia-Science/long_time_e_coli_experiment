library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


############

#df <- fread('geno_pred/plink/gwas_MAC250_bslmm_probit_top10.ld', header = T) %>% arrange(-R2)
#top_hits_file <- 'tables/gwas_bslmm_probit_top10_markers.txt'


#snakemake file args
ld_data_file <- snakemake@input[["ld_data"]]
top_hits_file <- snakemake@input[["top_hits"]]

output_heatmap_fig_file <- snakemake@output[["ld_heatmap_figure"]]


#########
#read in ld data
df <- fread(ld_data_file, header = T) %>% arrange(-R2)

#read in data on marker effects and phenotypes they are related to, focus only on ampicillin markers
top_hits <- fread(top_hits_file, header = T) %>%
    select(rs, rank, phenotype) %>% filter(phenotype == 'ampicillin') %>%
    mutate(marker_name = paste("Marker", rank)) %>%
    select(-rank)




#connect LD and marker data
df2 <- df %>% right_join(.,top_hits, by = c('SNP_A' = 'rs')) %>%
    rename(phenochr1 = phenotype) %>%
        right_join(.,top_hits, by = c('SNP_B' = 'rs')) %>%
            rename(phenochr2 = phenotype)




#select columns needed to make LD heatmap among markers
dfheat_single <- df2 %>% select(SNP_A, SNP_B, R2, marker_name.x, marker_name.y) %>%
        rename(Marker_A = marker_name.x, Marker_B = marker_name.y)

#same as above but flipping the identities of Marker A and Marker B to create the mirror image of matrix values for heatmap
dfheat_mirror <- df2 %>% select(SNP_A, SNP_B, R2, marker_name.x, marker_name.y) %>%
    rename(Marker_A = marker_name.y, Marker_B = marker_name.x) %>%
    rename(SNP_B_temp = SNP_A) %>% rename(SNP_A = SNP_B) %>% rename(SNP_B = SNP_B_temp)

#combine mirror and regular set of LD values for heatmap
dfheat <- rbind(dfheat_single, dfheat_mirror) %>% filter(!is.na(R2))


#convert marker identities to ordered factor so they are ordered from largest to smallest effect size
dfheat$Marker_A_ordered <- factor(dfheat$Marker_A, ordered=TRUE, levels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4",
                                                                            "Marker 5", "Marker 6", "Marker 7", "Marker 8",
                                                                            "Marker 9", "Marker 10"))
dfheat$Marker_B_ordered <- factor(dfheat$Marker_B, ordered=TRUE, levels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4",
                                                                            "Marker 5", "Marker 6", "Marker 7", "Marker 8",
                                                                            "Marker 9", "Marker 10"))

#plot LD heatmap among top10 ampicillin markers
plheat <- ggplot(data = dfheat, aes(x=Marker_A_ordered, y=Marker_B_ordered, fill=R2)) +
        geom_tile() + scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        scale_fill_gradientn(limits = c(0,0.75), colours=c("#341E60", "#A96789", "#F5DFB2"))+
        xlab("First marker") + ylab("Second marker") + guides(fill=guide_legend(title=expression(LD(~r^2))))

#save plot
#ggsave('figs/Fig6_LDheatmap_500.png', plheat, width = 6, height = 5)
ggsave(output_heatmap_fig_file, plheat, width = 6, height = 5)
