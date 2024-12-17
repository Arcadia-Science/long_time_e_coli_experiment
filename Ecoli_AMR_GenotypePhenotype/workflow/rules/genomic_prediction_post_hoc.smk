

#combines genomic pred output across phenotypes into one file, and creates a phenotype column from stripped file name
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


#plot top 10 largest marker effects
#generate site lists for bcftools and plink filtering of top10 markers by phenotype
#output table of top10 markers, and their gene annotations from pangenome (if available)
rule plot_gwas_results:
    conda: '../envs/popgenR.yaml'
    input:
        gwas_results_combined = expand("../geno_pred/output/gwas_MAC{mac}_bslmm_probit_combined.param.txt", mac = MAC),
        gene_names = rules.prep_gene_names.output.gene_names
    output:
        gwas_hits_fig = '../figs/Fig4_gwas_tophits_effects.png',
        gwas_top10_table = '../tables/gwas_bslmm_probit_top10_markers.txt',
        gwas_top10_plink_sites = '../geno_pred/plink/gwas_probit_top10_markers_ranges.txt',
        gwas_top10_sites = '../geno_pred/plink/gwas_probit_top10_markers_sites.txt'
    script:
        "scripts/gwas_results_presence_absence.R"


#subsequent analyses depend on specific markers uncovered during genomic prediction (gyrA248, gyrA259, parC239)
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


#generate gene tree for gyrA
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


#get genotypes for all SNPs in top10 marker lists from genomic prediction
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


#plot ciprofloxacin post hoc analyses results
#plot gene tree for gyrA locus and add ciprofloxacin resistance phenotype, and marker allele for markers gyrA248, gyrA259, parC239 (Fig5a)
#plot ciprofloxacin resistance level by genotype at 3 markers (Fig5b)
#run logistic regression of 3 markers vs ciprlfoxacin resistance
rule plot_ciprofloxacin_posthoc_results:
    conda: '../envs/popgenR.yaml'
    input:
        gyrA_tree = rules.generate_gyrA_tree.output.gyrA_tree,
        metadata_formatted = rules.format_metadata.output.metadata_formatted,
        marker_genotypes = rules.filter_gwas_top10_marker_genotypes.output.gts_gwas_top10_markers
    output:
        gyrA_tree_fig = '../figs/Fig5a_gyrA_tree_markers.png',
        cipro_marker_boxplot = '../figs/Fig5b_Boxplot_ciprofloxacin_marker_effects.png',
        cipro_marker_logistic_regression = '../tables/logistic_regression_ciprofloxacin.txt'
    script:
        "scripts/ciprofloxacin_posthoc_analyses.R"


#Time calibrate species tree for downstreaam corHMM analysis
rule time_calibrate_species_tree:
    conda: '../envs/popgenR.yaml'
    input:
        alignment_genome = rules.generate_wg_alignment.output.alignment_genome,
        uncalibrated_tree = rules.generate_species_tree.output.tree_genome,
        time_calibration_data = rules.format_metadata.output.time_calibration_data
    output:
        calibrated_tree_genome = '../tree/alignment/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAC10_refedit_subsample.min4.phy.timetree.nwk'
    shell:
        """
        iqtree -s {input.alignment_genome} \
        --date {input.time_calibration_data} -te {input.uncalibrated_tree} -m GTR+ASC --redo-tree
        """


#run corHMM to estimate transition rates between allelic states of 3 key ciprlofoxacin markers gyrA248, gyrA259, parC239
rule run_corHMM_ciprofloxacin_markers:
    conda: '../envs/popgenR.yaml'
    input:
        gyrA_tree_time_calibrated = rules.time_calibrate_species_tree.output.calibrated_tree_genome,
        marker_genotypes = rules.filter_gwas_top10_marker_genotypes.output.gts_gwas_top10_markers
    output:
        corhmm_ciprofloxaxin_markers = '../tables/corHMM_ciprofloxacin_markers.RDS'
    script:
        "scripts/ciprofloxacin_corhmm.R"


#plink calculate LD among all top10 markers in dataset
#first extract top10 markers for all phenos from bed file generates during genomic prediction then run all vs all ld calculation
rule plink_calculate_ld_gwas_markers:
    conda: '../envs/popgenR.yaml'
    input:
        input_bed = rules.merge_plink_input.output.combined_bed,
        gwas_top10_plink_sites = rules.plot_gwas_results.output.gwas_top10_plink_sites,
    output:
        output_bed_tophits = "../geno_pred/plink/gwas_MAC{mac}_bslmm_probit_top10.bed",
        plink_ld = "../geno_pred/plink/gwas_MAC{mac}_bslmm_probit_top10.ld"
    shell:
        """
        plink --bfile ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence  --make-bed --allow-extra-chr \
        --extract range  {input.gwas_top10_plink_sites} \
        --out ../geno_pred/plink/gwas_MAC{wildcards.mac}_bslmm_probit_top10

        plink --bfile ../geno_pred/plink/gwas_MAC{wildcards.mac}_bslmm_probit_top10 --r2 --allow-extra-chr \
        --ld-window-r2 0 --inter-chr --out ../geno_pred/plink/gwas_MAC{wildcards.mac}_bslmm_probit_top10
        """


#generate LD heatmap specifically for ampicillin top10 markers (Fig6)
rule ld_heatmap_ampicillin_markers:
    conda: '../envs/popgenR.yaml'
    input:
        ld_data = expand("../geno_pred/plink/gwas_MAC{mac}_bslmm_probit_top10.ld", mac=MAC),
        top_hits = rules.plot_gwas_results.output.gwas_top10_table
    output:
        ld_heatmap_figure = '../figs/Fig6_ampicillin_top10markers_ld_heatmap.png'
    script:
        "scripts/ld_heatmap_genopred.R"
