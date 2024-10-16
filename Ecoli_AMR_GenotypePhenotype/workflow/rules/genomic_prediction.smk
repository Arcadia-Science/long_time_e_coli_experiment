


wildcard_constraints:
    mac="|".join([str(m) for m in MAC]),
    pheno="|".join(PHENOTYPES),



#filter variants
rule filter_variants_genopred:
    input:
        vcf = "vcf_files/annotated_output.vcf.gz",
        outlier_samps = rules.identify_outlier_samps.output.outlier_samps
    output:
        vcf_genopred = "geno_pred/annotated_output_MAF{mac}_remoutliers.vcf"
    conda:
        "bcf"
    shell:
        "bcftools view {input.vcf} -S {input.outlier_samps} | \
        bcftools filter -i 'AC > {wildcards.mac}' | \
        sed 's/\.\:\./0\:255\,0/g' > {output.vcf_genopred}"



#generate plink files for SNP dataset
rule plink_generate_snp:
    input:
        vcf_genopred = rules.filter_variants_genopred.output.vcf_genopred
    output:
        snp_fam = "geno_pred/plink/annotated_output_MAC{mac}.fam",
        snp_bed = "geno_pred/plink/annotated_output_MAC{mac}.bed",
        snp_bim = "geno_pred/plink/annotated_output_MAC{mac}.bim"
    conda:
        "popgenR"
    shell:
        "plink --vcf {input.vcf_genopred} --make-bed --out geno_pred/plink/annotated_output_MAC{wildcards.mac} --allow-extra-chr"


#prepare SNP plink files for merging to presence absence data
rule edit_plink_siteid_snp:
    input:
        snp_bim = rules.plink_generate_snp.output.snp_bim,
    output:
        snp_bim_edited = 'geno_pred/plink/annotated_output_MAC{mac}.bim'
    shell:
        """
        awk '$2=$1' {input.snp_bim} |awk '{gsub(/_.*/,"",$2);print}' |awk '$2=NR $2 $4' > geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim
        mv geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim {output.snp_bim_edited}
        """


#prepare presence absence plink files for merging, need unique ID's for all sites for plink merge
#first adds bogus position, then generates a unique ID based on contig names and 'presabs' string
#copy over bed and fam, they don't need editing but to keep new files together
rule edit_plink_siteid_presence_absence:
    input:
        presence_absence_bim = 'presence_absence_data/presence_absence_all.bim',
        presence_absence_map = 'presence_absence_data/presence_absence_all.map',
        presence_absence_bed = 'presence_absence_data/presence_absence_all.bed',
        presence_absence_fam = 'presence_absence_data/presence_absence_all.fam'
    output:
        presence_absence_bim_edited = 'geno_pred/plink/presence_absence_all_edited.bim',
        presence_absence_map_edited = 'geno_pred/plink/presence_absence_all_edited.map',
        presence_absence_bed_edited = 'geno_pred/plink/presence_absence_all_edited.bed',
        presence_absence_fam_edited = 'geno_pred/plink/presence_absence_all_edited.fam'
    shell:
        """
        awk '{{gsub(/.*/,"696969",$4);print}}' {input.presence_absence_bim}  | \
        awk '$2=$1 "presabs"' > {output.presence_absence_bim_edited}

        awk '{{gsub(/.*/,"69696969",$3);print}}' {input.presence_absence_map} | \
        awk '$2=$1 "presabs"' > {output.presence_absence_map_edited}

        cp {input.presence_absence_bed}  {output.presence_absence_bed_edited}
        cp {input.presence_absence_fam}  {output.presence_absence_fam_edited}
        """




rule merge_plink_input:
    input:
        snp_bim = rules.edit_plink_siteid_snp.output.snp_bim_edited,
        snp_bed = rules.plink_generate_snp.output.snp_bed,
        snp_fam = rules.plink_generate_snp.output.snp_fam,
        presence_absence_bim = rules.edit_plink_siteid_presence_absence.output.presence_absence_bim_edited,
        presence_absence_map = rules.edit_plink_siteid_presence_absence.output.presence_absence_map_edited,
        presence_absence_bed = rules.edit_plink_siteid_presence_absence.output.presence_absence_bed_edited,
        presence_absence_fam = rules.edit_plink_siteid_presence_absence.output.presence_absence_fam_edited
    output:
        combined_bed = 'geno_pred/plink/annotated_output_MAC{mac}_presence_absence.bed',
        combined_bim = 'geno_pred/plink/annotated_output_MAC{mac}_presence_absence.bim',
        combined_fam = 'geno_pred/plink/annotated_output_MAC{mac}_presence_absence.fam'
    shell:
        """
        plink --bfile geno_pred/plink/presence_absence_all_edited --bmerge geno_pred/plink/annotated_output_MAC{wildcards.mac} \
        --allow-extra-chr --allow-no-sex --make-bed --out geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence
        """




#generate phenotype files
rule edit_plink_fam:
    input:
        rules.merge_plink_input.output.combined_fam
    output:
        fam_pheno_edit = "geno_pred/plink/annotated_output_MAC{mac}_{pheno}_presence_absence.fam"
    params:
        "{pheno}"
    conda:
        "popgenR"
    script:
        "geno_pred/gwas_plink_fam_edit.R"



#generate GEMMA relatedness matrix
rule run_GEMMA_relmatrix:
    input:
        rules.merge_plink_input.output.combined_bed
    output:
        rel_matrix = "output/annotated_output_MAC{mac}_presence_absence.cXX.txt"
    conda:
        "popgenR"
    shell:
        """
        sed -i 's/-9/1/g' {input}
        gemma -bfile {input}  -gk 1 -o annotated_output_MAC{wildcards.mac}_presence_absence
        """




#create temp gemma input files
rule GEMMA_temp_files:
    input:
        input_bed = rules.merge_plink_input.output.combined_bed,
        input_bim = rules.merge_plink_input.output.combined_bim
    output:
        bed_pheno_edit =  "geno_pred/plink/annotated_output_MAC{mac}_{pheno}presence_absence.bed",
        bim_pheno_edit =  "geno_pred/plink/annotated_output_MAC{mac}_{pheno}presence_absence.bim"
    params:
        "{pheno}"
    conda:
        "popgenR"
    shell:
        """
        cp {input.input_bed} {output.bed_pheno_edit}
        cp {input.input_bim} {output.bim_pheno_edit}
        """

#run GEMMA genopred probit
rule run_GEMMA_gwas_probit:
    input:
        genopred_bed = rules.GEMMA_temp_files.output.bed_pheno_edit,
        genopred_bim = rules.GEMMA_temp_files.output.bim_pheno_edit,
        genopred_fam = rules.edit_plink_fam.output.fam_pheno_edit,
        rel_matrix = rules.run_GEMMA_relmatrix.output.rel_matrix
    output:
        "geno_pred/output/gwas_MAC{mac}_{pheno}_bslmm_probit.param.txt"
    params:
        "{pheno}"
    conda:
        "popgenR"
    threads:10
    shell:
        """
        gemma -bfile geno_pred/plink/annotated_output_MAC{wildcards.mac}_{wildcards.pheno}presence_absence. -bslmm 3  -o gwas_MAC{wildcards.mac}_{wildcards.pheno}_bslmm_probit -k {input.rel_matrix}
        """
