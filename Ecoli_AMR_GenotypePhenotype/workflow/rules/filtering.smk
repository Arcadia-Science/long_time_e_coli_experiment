
#identify contigs shared by all ECOR strains, will be used to filter downstream
rule get_good_contigs:
    input:
        'pangenome/whole_pan_ecor_presence_absence.csv'
    output:
        'pangenome/ecor_shared_contigs.txt'
    script:
        "../../popgen_scripts/pan_genome.R"


rule filter_biallelic_annotated_snps:
    input:
        vcf = 'vcf_files/annotated_output.vcf.gz',
        contigs = 'pangenome/ecor_shared_contigs.txt'
    output:
        'vcf_files/annotated_output_biallelic_goodcontigs.vcf.gz'
    shell:
        """
        cat <(zgrep -e '#'{input.vcf}) \
        <(zgrep -Ev '#' {input.vcf} | zgrep -f {input.contigs} ) | zgrep -E 'synonymous|missense|HIGH|#' | \
         bcftools filter -e "TYPE = "indel""  | bcftools view -M2 -m2 -v snps -Oz  > {output}
        """
