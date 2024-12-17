

#identify contigs shared by all ECOR strains, will be used to filter downstream
rule get_good_contigs:
    conda: '../envs/popgenR.yaml'
    input:
        '../pangenome/whole_pan_ecor_presence_absence.csv'
    output:
        '../pangenome/ecor_shared_contigs.txt'
    script:
        "scripts/pan_genome.R"


#filter to onyl biallelic SNPs of 3 common classes on good contigs
rule filter_biallelic_annotated_snps:
    conda: '../envs/bcf.yaml'
    input:
        vcf = '../vcf_files/annotated_output.vcf.gz',
        contigs = '../pangenome/ecor_shared_contigs.txt'
    output:
        vcf_bi = '../vcf_files/annotated_output_biallelic_goodcontigs.vcf.gz'
    shell:
        """
        cat <(zgrep -e '#' {input.vcf}) \
        <(zgrep -Ev '#' {input.vcf} | zgrep -f {input.contigs} ) | zgrep -E 'synonymous|missense|HIGH|#' | \
        bcftools filter -e 'TYPE = "indel"'  | bcftools view -M2 -m2 -v snps -Oz  > {output.vcf_bi}
        """


#extract quadrupleton genotypes for synonymous SNPs to identify outlier samples in dataset
rule quads_analysis_input:
    conda: '../envs/bcf.yaml'
    input:
        vcf = rules.filter_biallelic_annotated_snps.output.vcf_bi
    output:
        gts_bi_syn = '../vcf_files/annotated_output_biallelic_goodcontigs_synonymous_quads_geno.txt'
    shell:
        """
        zgrep -E 'synonymous|#' {input.vcf} | \
        bcftools filter -i 'INFO/AC = 4'| \
        bcftools query -H -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t[\t%GT]\n" -o {output.gts_bi_syn}
        """


rule identify_outlier_samps:
    conda: '../envs/popgenR.yaml'
    input:
        vcf = rules.quads_analysis_input.output.gts_bi_syn
    output:
        outlier_samps = '../pangenome/outlier_samples.txt'
    script:
        "scripts/quadrupleton_stats.R"


rule filter_outlier_samples:
    conda: '../envs/bcf.yaml'
    input:
        vcf_bi = rules.filter_biallelic_annotated_snps.output.vcf_bi,
        outlier_samps = rules.identify_outlier_samps.output.outlier_samps
    output:
        vcf_outliers_removed = '../vcf_files/annotated_output_biallelic_goodcontigs_remoutliers.vcf.gz'
    shell:
        """
        bcftools view -S ^{input.outlier_samps} {input.vcf_bi}  -Oz  > {output.vcf_outliers_removed}
        """


rule generate_sfs_input:
    conda: '../envs/bcf.yaml'
    input:
        vcf_bi = rules.filter_biallelic_annotated_snps.output.vcf_bi,
        vcf_remoutlier = rules.filter_outlier_samples.output.vcf_outliers_removed
    output:
        sfs_allsamps = '../vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts.txt',
        sfs_outliers_removed = '../vcf_files/annotated_output_biallelic_goodcontigs_SFS_counts_remoutliers.txt'
    shell:
        """
        bcftools query  {input.vcf_bi} -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t%ANN\n" -o {output.sfs_allsamps}
        bcftools query  {input.vcf_remoutlier} -f "%CHROM\t%POS\t%ALT\t%AC\t%DP\t%ANN\n" -o {output.sfs_outliers_removed}
        """


rule generate_synonymous_filtered:
    conda: '../envs/bcf.yaml'
    input:
        vcf_remoutlier = rules.filter_outlier_samples.output.vcf_outliers_removed
    output:
        vcf_syn_popgen = '../vcf_files/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAC10_refedit.vcf'
    shell:
        """
         zgrep -E 'synonymous|#' {input.vcf_remoutlier} | \
         bcftools filter -i 'AC > 9'| sed 's/\.\:\./0\:255\,0/g'  > {output.vcf_syn_popgen}
        """


rule generate_popgen_input:
    conda: '../envs/bcf.yaml'
    input:
        vcf_syn_popgen = rules.generate_synonymous_filtered.output.vcf_syn_popgen
    output:
        vcf_syn_popgen_subsampled = '../tree/alignment/annotated_output_biallelic_goodcontigs_remoutliers_synonymous_MAC10_refedit_subsample.vcf'
    shell:
        """
        cat <(grep -e '#' {input.vcf_syn_popgen} ) <(grep -Ev '#'  {input.vcf_syn_popgen}| sed -n '0~10p') > {output.vcf_syn_popgen_subsampled}
        """
