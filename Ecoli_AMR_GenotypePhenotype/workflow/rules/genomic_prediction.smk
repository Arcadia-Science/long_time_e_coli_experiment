

wildcard_constraints:
    mac="|".join([str(m) for m in MAC]),
    pheno="|".join(PHENOTYPES),



#filter variants based on MAC cutoff defined
rule filter_variants_genopred:
    conda: '../envs/bcf.yaml'
    input:
        vcf = "../vcf_files/annotated_output.vcf.gz",
        outlier_samps = rules.identify_outlier_samps.output.outlier_samps
    output:
        vcf_genopred = "../geno_pred/annotated_output_MAF{mac}_remoutliers.vcf"
    shell:
        """
        bcftools view {input.vcf} -S ^{input.outlier_samps} | \
        bcftools filter -i 'AC > {wildcards.mac}' | \
        sed 's/\.\:\./0\:255\,0/g' > {output.vcf_genopred}
        """



#generate plink files for MAC filtered SNP dataset
rule plink_generate_snp:
    conda: '../envs/popgenR.yaml'
    input:
        vcf_genopred = rules.filter_variants_genopred.output.vcf_genopred
    output:
        snp_fam = "../geno_pred/plink/annotated_output_MAC{mac}.fam",
        snp_bed = "../geno_pred/plink/annotated_output_MAC{mac}.bed",
        snp_bim = "../geno_pred/plink/annotated_output_MAC{mac}_unedited.bim"
    shell:
        """
        plink --vcf {input.vcf_genopred} --make-bed --out ../geno_pred/plink/annotated_output_MAC{wildcards.mac} --allow-extra-chr
        mv ../geno_pred/plink/annotated_output_MAC{wildcards.mac}.bim ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_unedited.bim
        """


#prepare SNP plink files for merging to presence absence data by adding a marker ID
rule edit_plink_siteid_snp:
    conda: '../envs/popgenR.yaml'
    input:
        snp_bim = rules.plink_generate_snp.output.snp_bim,
    output:
        snp_bim_edited = '../geno_pred/plink/annotated_output_MAC{mac}.bim'
    shell:
        """
        awk '$2=$1' {input.snp_bim} |awk '{{gsub(/_.*/,"",$2);print}}' |awk '$2=NR $2 $4' > ../geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim
        mv ../geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim {output.snp_bim_edited}
        """


#prepare presence absence plink files for merging, need unique ID's for all sites for plink merge
#first adds bogus position (69696969), then generates a unique ID based on contig names and 'presabs' string
#copy over bed and fam, they don't need editing but to keep new files together
rule edit_plink_siteid_presence_absence:
    conda: '../envs/popgenR.yaml'
    input:
        presence_absence_bim = '../presence_absence_data/presence_absence_all.bim',
        presence_absence_map = '../presence_absence_data/presence_absence_all.map',
        presence_absence_bed = '../presence_absence_data/presence_absence_all.bed',
        presence_absence_fam = '../presence_absence_data/presence_absence_all.fam'
    output:
        presence_absence_bim_edited = '../geno_pred/plink/presence_absence_all_edited.bim',
        presence_absence_map_edited = '../geno_pred/plink/presence_absence_all_edited.map',
        presence_absence_bed_edited = '../geno_pred/plink/presence_absence_all_edited.bed',
        presence_absence_fam_edited = '../geno_pred/plink/presence_absence_all_edited.fam'
    shell:
        """
        awk '{{gsub(/.*/,"696969",$4);print}}' {input.presence_absence_bim}  | \
        awk '$2=$1 "presabs"' > {output.presence_absence_bim_edited}

        awk '{{gsub(/.*/,"69696969",$3);print}}' {input.presence_absence_map} | \
        awk '$2=$1 "presabs"' > {output.presence_absence_map_edited}

        cp {input.presence_absence_bed}  {output.presence_absence_bed_edited}
        cp {input.presence_absence_fam}  {output.presence_absence_fam_edited}
        """



#merge the SNP bed and presence/absence marker bed files
rule merge_plink_input:
    conda: '../envs/popgenR.yaml'
    input:
        snp_bim2 = rules.edit_plink_siteid_snp.output.snp_bim_edited,
        snp_bed = rules.plink_generate_snp.output.snp_bed,
        snp_fam = rules.plink_generate_snp.output.snp_fam,
        presence_absence_bim = rules.edit_plink_siteid_presence_absence.output.presence_absence_bim_edited,
        presence_absence_map = rules.edit_plink_siteid_presence_absence.output.presence_absence_map_edited,
        presence_absence_bed = rules.edit_plink_siteid_presence_absence.output.presence_absence_bed_edited,
        presence_absence_fam = rules.edit_plink_siteid_presence_absence.output.presence_absence_fam_edited
    output:
        combined_bed = '../geno_pred/plink/annotated_output_MAC{mac}_presence_absence.bed',
        combined_bim = '../geno_pred/plink/annotated_output_MAC{mac}_presence_absence.bim',
        combined_fam = '../geno_pred/plink/annotated_output_MAC{mac}_presence_absence.fam'
    shell:
        """
        plink --bfile ../geno_pred/plink/presence_absence_all_edited --bmerge ../geno_pred/plink/annotated_output_MAC{wildcards.mac} \
        --allow-extra-chr --allow-no-sex --make-bed --out ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence
        """




#generate phenotype plink files (.fam) by adding susceptible(0) and resistant/intermediate(1) column from metadata file
rule edit_plink_fam:
    conda: '../envs/popgenR.yaml'
    input:
        rules.merge_plink_input.output.combined_fam,
        rules.format_metadata.output.metadata_formatted
    output:
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.fam"
    params:
        "{pheno}"
    conda:
        "popgenR"
    script:
        "scripts/gwas_plink_fam_edit.R"



#generate GEMMA relatedness matrix to save time since same matrix will be used for all phenotypes
rule run_GEMMA_relmatrix:
    conda: '../envs/popgenR.yaml'
    input:
        rules.merge_plink_input.output.combined_bed
    output:
        "../geno_pred/output/annotated_output_MAC{mac}.cXX.txt"
    params:
        genopred_directory = '../geno_pred/'
    shell:
        """
        cd {params.genopred_directory}
        sed -i 's/-9/1/g' ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence.fam
        gemma -bfile ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence  -gk 1 -o annotated_output_MAC{wildcards.mac}
        """




#create temporary phenotype specific (non fam) gemma input files
rule GEMMA_temp_files:
    conda: '../envs/popgenR.yaml'
    input:
        rules.merge_plink_input.output.combined_bed,
        rules.merge_plink_input.output.combined_bim
    output:
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.bed",
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.bim"
    params:
        "{pheno}"
    shell:
        """
        cp ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence.bed ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence_{wildcards.pheno}.bed
        cp ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence.bim ../geno_pred/plink/annotated_output_MAC{wildcards.mac}_presence_absence_{wildcards.pheno}.bim
        """

#run genomic prediction GEMMA genopred probit model
rule run_GEMMA_gwas_probit:
    conda: '../envs/popgenR.yaml'
    input:
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.bed",
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.bim",
        "../geno_pred/plink/annotated_output_MAC{mac}_presence_absence_{pheno}.fam",
        "../geno_pred/output/annotated_output_MAC{mac}.cXX.txt"
    output:
        "../geno_pred/output/gwas_MAC{mac}_{pheno}_bslmm_probit.param.txt"
    params:
        "{pheno}",
        genopred_directory = '../geno_pred/'
    threads:13
    resources:
        mem_mb=10000
    shell:
        """
        cd {params.genopred_directory}
        gemma -bfile plink/annotated_output_MAC{wildcards.mac}_presence_absence_{wildcards.pheno} -bslmm 3  -o gwas_MAC{wildcards.mac}_{wildcards.pheno}_bslmm_probit -k output/annotated_output_MAC{wildcards.mac}.cXX.txt
        """
