

#combines geno pred output across phenotypes into one file, and creates a phenotype column from stripped file name
rule GEMMA_combine_output:
    conda: '../envs/popgenR.yaml'
    input:
        expand("../geno_pred/output/gwas_MAC{mac}_{pheno}_bslmm_probit.param.txt", mac = MAC, pheno = PHENOTYPES)
    output:
        gwas_results_combined = "../geno_pred/output/gwas_MAC{mac}_bslmm_probit_combined.param.txt"
    shell:
        """
        awk '{{print $0, FILENAME}}' {input} | \
        awk '{{gsub(/\.\.\/geno.*MAC[0-9]*_|_bslmm.*/,"",$NF);print}}'  | sed '/^chr.*$/{{x;/^$/!d;g;}}'  >> {output}
        """

#scrape gene names from pangenome reference fasta to help interpret genomic prediction results
rule prep_gene_names:
    conda: '../envs/popgenR.yaml'
    input:
        '../pangenome/whole_pangenome.fasta'
    output:
        gene_names = "../pangenome/gene_names.txt"
    shell:
        """
        grep -e '>' {input} | sed 's/>//g' > {output}
        """

#plot top 10 largest marker effects and generate marker list for further plink analyses
rule plot_gwas_results:
    conda: '../envs/popgenR.yaml'
    input:
        gwas_results_combined = expand("../geno_pred/output/gwas_MAC{mac}_bslmm_probit_combined.param.txt", mac = MAC),
        gene_names = rules.prep_gene_names.output.gene_names
    output:
        gwas_hits_fig = '../figs/gwas_tophits_effects.png',
        gwas_top10_table = '../tables/gwas_bslmm_probit_top10_markers.txt',
        gwas_top10_plink_sites = '../geno_pred/plink/gwas_probit_top10_markers_ranges.txt',
        gwas_top10_sites = '../geno_pred/plink/gwas_probit_top10_markers_sites.txt'
    script:
        "scripts/gwas_results_presence_absence.R"

#get gyrA data for gene tree building
rule generate_gyrA_tree_input:
    conda: '../envs/bcf.yaml'
    input:
        input_vcf = rules.filter_outlier_samples.output.vcf_outliers_removed
    output:
        gyrA_vcf = '../tree/alignment/gyrA_remoutliers.vcf',
        gyrA_vcf_edited = '../tree/alignment/gyrA_remoutliers_edited.vcf',
        gyrA_alignment = '../tree/alignment/gyrA_remoutliers_edited.min4.phy'
    shell:
        """
        rules/scripts/generate_gyrA_tree_input.sh {input.input_vcf}  {output.gyrA_vcf} {output.gyrA_vcf_edited}
        """


rule generate_gyrA_tree:
    conda: '../envs/popgenR.yaml'
    input:
        alignment_genome = rules.generate_gyrA_tree_input.output.gyrA_alignment
    output:
        gyrA_tree = '../tree/alignment/gyrA_remoutliers_edited.min4.phy.treefile',
    shell:
        """
        iqtree -s {input.alignment_genome}  -m GTR+ASC
        """



rule filter_gwas_top10_marker_genotypes:
    conda: '../envs/bcf.yaml'
    input:
        vcf = '../vcf_files/annotated_output.vcf.gz',
        marker_sites = rules.plot_gwas_results.output.gwas_top10_sites,
        outlier_samps = rules.identify_outlier_samps.output.outlier_samps
    output:
        gts_gwas_top10_markers = '../geno_pred/gwas_top10_marker_genotypes.txt'
    shell:
        """
        bcftools view -T {input.marker_sites} -S ^{input.outlier_samps} {input.vcf} | \
        bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP[\t%GT]\n" -o {output.gts_gwas_top10_markers}
        """



rule plot_ciprofloxacin_posthoc_results:
    conda: '../envs/popgenR.yaml'
    input:
        gyrA_tree = rules.generate_gyrA_tree.output.gyrA_tree,
        metadata_formatted = rules.format_metadata.output.metadata_formatted,
        marker_genotypes = rules.filter_gwas_top10_marker_genotypes.output.gts_gwas_top10_markers
    output:
        gyrA_tree_fig = '../figs/gyrA_tree_markers.png',
        cipro_marker_boxplot = '../figs/Boxplot_ciprofloxacin_marker_effects.png',
        cipro_marker_logistic_regression = '../tables/logistic_regression_ciprofloxacin.txt'
    script:
        "scripts/ciprofloxacin_posthoc_analyses.R"
