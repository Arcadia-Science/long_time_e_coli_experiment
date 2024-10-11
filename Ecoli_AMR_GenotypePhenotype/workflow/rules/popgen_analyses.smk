

#generate combined metdata files for downstream analyses
rule format_metadata:
    input:
        SRA_GCA = 'pangenome/pangenome_genomes_SRA_GCA.csv',
        phenotype_matrix = 'phenotype_matrix/phenotype_matrix_08302024.csv',
        SRA_genome_name = 'phenotype_matrix/SRA_to_genome_name.csv',
        strains_metadata = 'strains_metadata.csv',
        outlier_samples = rules.identify_outlier_samps.output.outlier_samps
    output:
        metadata_formatted = 'strains_metadata_phenotypes_full.txt',
        sample_list = 'strains_sample_list.txt',
        time_calibration_data = 'strains_year_collection.txt',
    script:
        "../../popgen_scripts/format_metadata.R"



#run PCA, generate PCA plots and multinomial models of PCA vs metadata features like country/year of isolation etc.
rule genotypic_PCA:
    input:
        popgen_syn_vcf = rules.generate_popgen_input.output.vcf_syn_popgen_subsampled,
        metadata_formatted = rules.format_metadata.output.metadata_formatted
    output:
        pca_fig = 'figs/genomic_PCA.txt',
        pca_model_output = 'tables/PCA_multinomial_model_fits.txt'
    script:
        "../../popgen_scripts/PCA.R"


#plot SFS
rule generate_SFS_plots:
    input:
        sfs_no_outlier_input = rules.generate_sfs_input.output.sfs_outliers_removed
    output:
        sfs_fig = 'figs/SFS_remoutliers.txt'
    script:
        "../../popgen_scripts/SFS.R"





#plot presence absence

#generate species tree and generate tree plots
