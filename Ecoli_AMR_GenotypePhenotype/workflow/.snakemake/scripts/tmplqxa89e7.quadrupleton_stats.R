
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('../vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt', "vcf" = '../vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt'),
    output = list('../pangenome/outlier_samples.txt', "outlier_samps" = '../pangenome/outlier_samples.txt'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("samples" = '../resources/glue_pc_sampleSheet.txt', "chromosomes" = '../resources/chromosome_file.txt', "gff" = '../resources/TrR.v5.renamed_reformated.gtf.gz', "vcf_prefix" = '/vcf_files', "pangenome_prefix" = '/pangenome'),
    rule = 'identify_outlier_samps',
    bench_iteration = as.numeric(NA),
    scriptdir = '/home/ubuntu/long_time_e_coli_experiment/Ecoli_AMR_GenotypePhenotype/workflow/rules/scripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#df <- fread('../vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt')
input_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]


df <- fread(input_file)

#funs to cleat bcftools column names
colClean_bracket <- function(x){ colnames(x) <- gsub(".*\\]", "", colnames(x)); x }
colClean_GT     <- function(x){ colnames(x) <- gsub("\\:.*", "", colnames(x)); x }
convert_numeric <- function(x){ as.numeric(x)}

df2 <- df %>% colClean_bracket(.) %>% colClean_GT(.) %>%
                mutate(snp_id = paste(CHROM, POS, sep = '_')) %>%
                relocate(snp_id, .before = CHROM) %>%
                select(-CHROM, -POS, -ALT, -AC, -DP, -V6) %>%
                mutate(across(!snp_id, convert_numeric))
head(df2)
#viz results
df3 <- df2 %>% select(-snp_id) %>% colSums(., na.rm=T)

#pull out outliers
outlier_samples <- df3 %>% as.data.frame(.) %>% filter(`.` > 500) %>% tibble:::rownames_to_column(.) %>% select(rowname)

write.table(outliers_samples, out_file, quote = F, col.names = F, row.names = T, sep = '\t' )

#write.table(outlier_samples, 'pangenome/outlier_samples.txt', quote = F, col.names = F, row.names = F, sep = '\t' )