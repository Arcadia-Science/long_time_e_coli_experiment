

#generate combined metdata files for downstream analyses
rule format_metadata:
    conda: '../envs/popgenR.yaml'
    input:
        SRA_GCA = '../pangenome/pangenome_genomes_SRA_GCA.csv',
        phenotype_matrix = '../phenotype_matrix/phenotype_matrix_08302024.csv',
        SRA_genome_name = '../phenotype_matrix/SRA_to_genome_name.csv',
        strains_metadata = '../strains_metadata.csv',
        outlier_samples = rules.identify_outlier_samps.output.outlier_samps
    output:
        metadata_formatted = '../strains_metadata_phenotypes_full.txt',
        sample_list = '../strains_sample_list.txt',
        time_calibration_data = '../strains_year_collection.txt',
    script:
        "scripts/format_metadata.R"



#run PCA, generate PCA plots and multinomial models of PCA vs metadata features like country/year of isolation etc.
rule genotypic_PCA:
    conda: '../envs/popgenR.yaml'
    input:
        popgen_syn_vcf = rules.generate_popgen_input.output.vcf_syn_popgen_subsampled,
        metadata_formatted = rules.format_metadata.output.metadata_formatted
    output:
        pca_fig = '../figs/genomic_PCA.png',
        pca_model_output = '../tables/PCA_multinomial_model_fits.txt'
    script:
        "scripts/PCA.R"


#plot SFS
rule generate_sfs_plots:
    conda: '../envs/popgenR.yaml'
    input:
        sfs_no_outlier_input = rules.generate_sfs_input.output.sfs_outliers_removed,
        presence_absence_input = '../presence_absence_data/presence_absence_sfs.txt',
        presence_absence_loci = '../presence_absence_data/presence_absence_all.bim'
    output:
        sfs_fig = '../figs/SFS_combined_remoutliers.png'
    script:
        "scripts/SFS.R"


#download repo for converting VCF to alignment file for phylogenomic scripts
rule clone_vcf2phylip:
    output:
        directory('../vcf2phylip'),
        '../vcf2phylip/vcf2phylip.py'
    shell:
        """
        git clone https://github.com/edgardomortiz/vcf2phylip.git ../vcf2phylip
        """


#generate alignment from vcf
rule generate_wg_alignment:
    input:
        '../vcf2phylip/vcf2phylip.py',
        popgen_vcf = rules.generate_popgen_input.output.vcf_syn_popgen_subsampled,
    output:
        alignment_genome = '../tree/alignment/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAC10_refedit_subsample.min4.phy',
    shell:
        """
        ../vcf2phylip/vcf2phylip.py --output-folder ../tree/alignment/ -i {input.popgen_vcf}
        """

#generate species tree
rule generate_species_tree:
    conda: '../envs/popgenR.yaml'
    input:
        alignment_genome = rules.generate_wg_alignment.output.alignment_genome
    output:
        tree_genome = '../tree/alignment/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAC10_refedit_subsample.min4.phy.treefile',
    shell:
        """
        iqtree -s {input.alignment_genome}  -m GTR+ASC
        """


#generate tree plots
rule generate_genome_phylogeny_plot:
    conda: '../envs/popgenR.yaml'
    input:
        tree_genome = rules.generate_species_tree.output.tree_genome,
        metadata_formatted = rules.format_metadata.output.metadata_formatted
    output:
        tree_genome_plot = '../figs/genome_tree_synonymous_MAC10.png',
    script:
        """
        "scripts/plot_genome_tree.R"
        """