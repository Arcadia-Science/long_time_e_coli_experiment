


wildcard_constraints:
    mac="|".join([str(m) for m in MAC]),
    pheno="|".join(PHENOTYPES),



#filter variants
rule filter_variants_genopred:
    conda: '../envs/bcf.yaml'
    input:
        vcf = "../vcf_files/annotated_output.vcf.gz",
        outlier_samps = rules.identify_outlier_samps.output.outlier_samps
    output:
        vcf_genopred = "../geno_pred/annotated_output_MAF{mac}_remoutliers.vcf"
    conda:
        "bcf"
    shell:
        "bcftools view {input.vcf} -S {input.outlier_samps} | \
        bcftools filter -i 'AC > {wildcards.mac}' | \
        sed 's/\.\:\./0\:255\,0/g' > {output.vcf_genopred}"



#generate plink files for SNP dataset
rule plink_generate_snp:
    conda: '../envs/popgenR.yaml'
    input:
        vcf_genopred = rules.filter_variants_genopred.output.vcf_genopred
    output:
        snp_fam = "../geno_pred/plink/annotated_output_MAC{mac}.fam",
        snp_bed = "../geno_pred/plink/annotated_output_MAC{mac}.bed",
        snp_bim = "../geno_pred/plink/annotated_output_MAC{mac}.bim"
    shell:
        "plink --vcf {input.vcf_genopred} --make-bed --out ../geno_pred/plink/annotated_output_MAC{wildcards.mac} --allow-extra-chr"


#prepare SNP plink files for merging to presence absence data
rule edit_plink_siteid_snp:
    conda: '../envs/popgenR.yaml'
    input:
        snp_bim = rules.plink_generate_snp.output.snp_bim,
    output:
        snp_bim_edited = '../geno_pred/plink/annotated_output_MAC{mac}.bim'
    shell:
        """
        awk '$2=$1' {input.snp_bim} |awk '{gsub(/_.*/,"",$2);print}' |awk '$2=NR $2 $4' > ../geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim
        mv ../geno_pred/plink/tmp_bim_snp_edited_{wildcards.mac}.bim {output.snp_bim_edited}
        """


#prepare presence absence plink files for merging, need unique ID's for all sites for plink merge
#first adds bogus position, then generates a unique ID based on contig names and 'presabs' string
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
        awk '{gsub(/.*/,"696969",$4);print}' {input.presence_absence_bim}  | \
        awk '$2=$1 "presabs"' > {output.presence_absence_bim_edited}

        awk '{gsub(/.*/,"69696969",$3);print}' {input.presence_absence_map} | \
        awk '$2=$1 "presabs"' > {output.presence_absence_map_edited}

        cp {input.presence_absence_bed}  {output.presence_absence_bed_edited}
        cp {input.presence_absence_fam}  {output.presence_absence_fam_edited}
        """




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




#generate phenotype files
rule edit_plink_fam:
    conda: '../envs/popgenR.yaml'
    input:
        "../geno_pred/plink/annotated_output_MAF{mac}.fam"
    output:
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.fam"
    params:
        "{pheno}"
    conda:
        "popgenR"
    script:
        "scripts/gwas_plink_fam_edit.R"



#generate GEMMA relatedness matrix
rule run_GEMMA_relmatrix:
    conda: '../envs/popgenR.yaml'
    input:
        "../geno_pred/plink/annotated_output_MAF{mac}.bed"
    output:
        "../geno_pred/output/annotated_output_MAF{mac}.cXX.txt"
    param:
        genopred_directory = '../geno_pred/'
    shell:
        """
        cd {params.genopred_directory}
        sed -i 's/-9/1/g' plink/annotated_output_MAF{wildcards.mac}.fam
        gemma -bfile plink/annotated_output_MAF{wildcards.mac}  -gk 1 -o annotated_output_MAF{wildcards.mac}
        """




#create temp gemma input files
rule GEMMA_temp_files:
    input:
        "../geno_pred/plink/annotated_output_MAF{mac}.bed",
        "../geno_pred/plink/annotated_output_MAF{mac}.bim"
    output:
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.bed",
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.bim"
    params:
        "{pheno}"
    conda:
        "popgenR"
    shell:
        """
        cp ../geno_pred/plink/annotated_output_MAF{wildcards.mac}.bed ../geno_pred/plink/annotated_output_MAF{wildcards.mac}_{wildcards.pheno}.bed
        cp ../geno_pred/plink/annotated_output_MAF{wildcards.mac}.bim ../geno_pred/plink/annotated_output_MAF{wildcards.mac}_{wildcards.pheno}.bim
        """

#run GEMMA genopred probit
rule run_GEMMA_gwas_probit:
    input:
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.bed",
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.bim",
        "../geno_pred/plink/annotated_output_MAF{mac}_{pheno}.fam",
        "../geno_pred/output/annotated_output_MAF{mac}.cXX.txt"
    output:
        "../geno_pred/output/gwas_MAF{mac}_{pheno}_bslmm_probit.param.txt"
    params:
        "{pheno}"
    conda:
        "popgenR"
    threads:10
    shell:
        """
        gemma -bfile ../geno_pred/plink/annotated_output_MAF{wildcards.mac}_{wildcards.pheno} -bslmm 3  -o gwas_MAF{wildcards.mac}_{wildcards.pheno}_bslmm_probit -k output/annotated_output_MAF{wildcards.mac}.cXX.txt
        """
