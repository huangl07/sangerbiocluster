#!/usr/bin/env nextflow
params.sample_tumor_pair="sample_tumor_pair"
params.wgs_result_dir="wgs_result_dir"
params.hla="yes"
params.hlalist="HLA-A"
params.outdir="./"
params.tumoronly="no"
def helpMessage(){
	log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --outdir    <dir>   output dir
    --ref_name <string> GRCH38/GRCH37
    --ref   <file>  reference genome fasta file
    --gff   <file>  reference gff file
    --chrlist   <file>  chrlist for draw
    --bamlist   <file>  bamlist for analyze
    --samplelist   <file>  samplelist for analyze
	--sample_tumor_pair   <file>  sample_tumor_pair for analyze
    --grouplist   <file>  grouplist for analyze
    --vcf <file> wgs result
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
pair_info=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info2=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info3=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info4=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info5=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info6=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info7=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
pair_info8=Channel.from(file(params.sample_tumor_pair))
              .splitCsv(header:false,sep:'\t')
gnomAD = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz"
dbsnp = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/GCF_000001405.39.gz"
cosmic = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/CosmicCodingMuts.vcf.gz"
dbNSFP_hg38 = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/dbNSFP4.3a.hg38.txt.gz"
dbNSFP_hg19 = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/dbNSFP4.3a.hg19.txt.gz"
change_chr_name_dbsnp = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/predisposing_genes_database/change_chr_name_dbsnp.txt"
snpsift = "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/snpEff/exec/snpsift"

wgs=Channel.from(file(params.sample_tumor_pair))
mhc_I_list = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/MHC_data/MHC_allele.list"
mhc_II_list = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/MHC_data/MHCII_allele.list"
gtf = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/ref/ref.gtf"
ref = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/ref/ref.fa"
ref_protein = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/ref/ref.protein.fa"
anno_summary = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/ref/anno.summary"
mhc_I = Channel.from(file(mhc_I_list))
              .splitCsv(header:false,sep:'\t')
mhc_II = Channel.from(file(mhc_II_list))
              .splitCsv(header:false,sep:'\t')
gtf_ch = Channel.from(file(gtf))
ref_ch = Channel.from(file(ref))

process tcga_mutect{
	publishDir "${params.outdir}/00.basic_anlysis/" , pattern: "*"
	tag "tcga_mutect"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        val(wgs) from wgs
    output:
        file "*"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    gunzip -c ${params.wgs_result_dir}/04.snpIndel/finalvcf/pop.snpindel.final.vcf.gz > pop.snpindel.final.vcf
    bash ${baseDir}/bin/00.tcga_mutect.sh -v pop.snpindel.final.vcf
    """    
}

process generate_full_vcf{
	publishDir "${params.outdir}/05.hf_mutation_gene_anlysis/full_vcf/" , pattern: "*"
	tag "generate_full_vcf"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control) from pair_info
    output:
        tuple val(sample),val(control),path("pop.snpindel.final.vcf") into full_vcf
        tuple val(sample),val(control),path("pop.snpindel.final.vcf") into full_vcf2
        tuple val(sample),val(control),path("pop.snpindel.final.vcf") into full_vcf3
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    gunzip -c ${params.wgs_result_dir}/04.snpIndel/finalvcf/pop.snpindel.final.vcf.gz > pop.snpindel.final.vcf
    """    
}
process generate_sample_control_vcf{
	publishDir "${params.outdir}/05.hf_mutation_gene_anlysis/sample_control_vcf/" , pattern: "*"
	tag "generate_sample_control_vcf"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control),path(full_vcf) from full_vcf
    output:
        tuple val(sample),val(control),path("${sample}.vcf"),path("${control}.vcf")into ms_vcf
        tuple val(sample),val(control) into ms_vcf2
        tuple path("${sample}.vcf"),path("${control}.vcf") into ms_vcf3
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    bcftools view ${full_vcf} -s ${sample} > sample_test.vcf
    less sample_test.vcf |grep -v "0/0"| grep -v "\\./\\."> ${sample}.vcf
    bcftools view ${full_vcf} -s ${control} > control_test.vcf
    less control_test.vcf | grep -v "0/0" | grep -v "\\./\\."> ${control}.vcf
    """    
}
process mutation_spectrum_pattern_hg38{
        publishDir "${params.outdir}/05.hf_mutation_gene_anlysis/mutation_spectrum/" , pattern: "*"
        tag "mutation_spectrum_pattern_hg38"
        queue "SANGERDEV"
        executor "slurm"
        cpus 8
        memory "32G"
        input:
            tuple val(sample),val(control),path(sample_vcf),path(control_vcf) from ms_vcf
        output:
            file "*"
        script:
        """
        source  ~/app/bioinfo/dna/new.rc
        Rscript ${baseDir}/bin/05.hi_freq_mutation/05.mutation_spectrum.R --samplefile ${sample_vcf} --controlfile ${control_vcf} --samplename ${sample} --controlname ${control} --refinfo BSgenome.Hsapiens.UCSC.hg38 --out ./${sample}_${control}
        Rscript ${baseDir}/bin/05.hi_freq_mutation/05.mutation_pattern.R --samplefile ${sample_vcf} --controlfile ${control_vcf} --samplename ${sample} --controlname ${control} --refinfo BSgenome.Hsapiens.UCSC.hg38 --signatures ${baseDir}/bin/05.hi_freq_mutation/COSMIC_v3.3.1_SBS_GRCh38.txt --out ./${sample}_${control}
        """  
}

process generate_sample_list{
    publishDir "${params.outdir}/05.hf_mutation_gene_anlysis/mutation_spectrum_merge/" , pattern: "*"
    tag "generate_sample_list"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "64G"
    input:
        val(sample) from ms_vcf2.collect()
        path(vcf) from ms_vcf3.collect()
    output:
        file "*"
    script:
    def all_sample = sample.collect().join(",")
    def all_vcf = vcf.collect().join(",")
    """
    source  ~/app/bioinfo/dna/new.rc
    Rscript ${baseDir}/bin/05.hi_freq_mutation/05.mutation_all.R --samplename ${all_sample} --samplefile  ${all_vcf} --refinfo BSgenome.Hsapiens.UCSC.hg38  --signatures ${baseDir}/bin/05.hi_freq_mutation/COSMIC_v3.3.1_SBS_GRCh38.txt --out ./
    convert -density 300 ./all_mutation_spectrum.pdf -quality 100 ./all_mutation_spectrum.png
    """    
}

process generate_maf{
    publishDir "${params.outdir}/04.clone_anlysis/clone_analysis_maf_vcf/" , pattern: "*"
    tag "generatemaf"
    queue "SANGERDEV"
    executor "slurm"
    cpus 4
    memory "8G"
    input:
        tuple val(sample),val(control),path(full_vcf) from full_vcf2
    output:
        tuple val("${sample}_${control}"),path("${sample}_${control}.maf")into maf_ch
        path("${sample}_${control}.mure.maf")into maf_ch1
        file("${sample}.${control}.vcf")
        tuple val("${sample}_${control}"),file("${sample}_${control}.vcf") into pre_gene_vcf
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    bcftools view -s ${sample},${control} ${full_vcf} > ${sample}_${control}.vcf
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\tGT:DP:AD[\\t%GT:%DP:%AD]\\n' ${sample}_${control}.vcf -o ${sample}.${control}.vcf
    perl ${baseDir}/bin/snptomaf_for_clone.pl ${sample}.${control}.vcf 100 0.05
    perl ${baseDir}/bin/snptomaf_for_mure.pl ${sample}.${control}.vcf 100 0.05
    """
}
process plot_mure{
    publishDir "${params.outdir}/07.significantly_mutated/mutation_relation/" , pattern: "*"
    tag "plotmure"
    queue "SANGERDEV"
    executor "slurm"
    cpus 2
    memory "4G"
    input:
        path(mafs) from maf_ch1.collect()
    output:
        file("*")
    script:
    def all_maf = mafs.collect().join(" ")
    """
    source  ~/app/bioinfo/dna/new.rc
    cat ${baseDir}/bin/maf_header_mure.txt ${all_maf}> all.maf
    Rscript ${baseDir}/bin/panorama_of_genomic_mutations.R --maf all.maf --out ./
    Rscript ${baseDir}/bin/mutation_relation_modified.R --maf all.maf --out ./
    """
}

process plot_infergenity{
    publishDir "${params.outdir}/04.clone_anlysis/clone_analysis/" , pattern: "*"
    tag "plotinfergenity"
    queue "SANGERDEV"
    executor "slurm"
    cpus 4
    memory "8G"
    input:
        tuple val(tumor_name),path(tumor_maf) from maf_ch
    output:
        file("*")
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    Rscript ${baseDir}/bin/clone_analysis.R --maf ${tumor_name}.maf --out ${tumor_name}
    Rscript ${baseDir}/bin/oncodrive.R --maf ${tumor_name}.maf --out ./${tumor_name}
    """
}

process sequenza{
    publishDir "${params.outdir}/01.cnl/sequenza_raw_data/" , pattern: "*"
    tag "sequenza"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(sample),val(control) from pair_info2
    output:
        tuple val(sample),val(control),file("*.txt.gz"),file("${sample}_${control}_small.seqz") into sequenza
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    sequenza-utils gc_wiggle -w 50 --fasta ${ref} |gzip > hg38.gc50Base.txt.gz
    sequenza-utils bam2seqz -n ${params.wgs_result_dir}/01.fastq_qc/cram/${control}.mkdup.cram -t ${params.wgs_result_dir}/01.fastq_qc/cram/${sample}.mkdup.cram --fasta ${ref} -gc hg38.gc50Base.txt.gz -o ${sample}_${control}.seqz
    sequenza-utils seqz_binning -w 50 -s ${sample}_${control}.seqz -o ${sample}_${control}_small.seqz
    """
}
process sequenza_plot{
    publishDir "${params.outdir}/01.cnl/sequenza_plot/" , pattern: "*"
    tag "sequenza_plot"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(sample),val(control),file(hg38),file(small_seqz) from sequenza
    output:
        file "*"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    mkdir ${sample}_${control}
    cd ${sample}_${control}
    Rscript ${baseDir}/bin/sequenza_plot.R --seqgz ../${small_seqz} --outpath ./ --samplename ${sample}_${control}
    convert -density 300 ${sample}_${control}_CP_contours.pdf -quality 100 ${sample}_${control}_CP_contours.png
    python /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/aCNViewer/code/aCNViewer.py -f ./ -t ${sample}_${control} --refBuild hg38 -w 2000000 -b /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/aCNViewer_DATA/bin/ --fileType Sequenza
    """
}

process filter_anno_vcf {
    publishDir "${params.outdir}/02.pre_gene/filtered_anno_vcf/", pattern: "${sample_name}.filtered_pool.vcf.gz"
    tag "filter_anno_vcf"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), file(splited_anno_vcf) from pre_gene_vcf
    output:
        tuple val(sample_name), path("${sample_name}.filtered_pool.vcf.gz") into filtered_pool_ch
        tuple val(sample_name), path("${sample_name}.filtered_pool_new.vcf.gz"), path("${sample_name}.filtered_pool_new.vcf.gz.tbi") into filtered_anno_vcf_ch1, filtered_anno_vcf_ch2, filtered_anno_vcf_ch3
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        ${snpsift} filter "(FILTER='PASS')" ${splited_anno_vcf} | bcftools view -O z > filtered_by_filter.vcf.gz
        ${snpsift} filter "(GEN[ALL].DP>10)" filtered_by_filter.vcf.gz | bcftools view -O z > filtered_by_dp.vcf.gz
        ${snpsift} filter "(ANN[*].IMPACT='HIGH')|(ANN[*].IMPACT='MODERATE')" filtered_by_dp.vcf.gz | bcftools view -O z > filtered_by_impact.vcf.gz
        bcftools filter -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY filtered_by_impact.vcf.gz -Oz -o ${sample_name}.filtered_pool.vcf.gz
        bcftools annotate -x INFO ${sample_name}.filtered_pool.vcf.gz -Oz -o ${sample_name}.filtered_pool_new.vcf.gz
        tabix -p vcf ${sample_name}.filtered_pool_new.vcf.gz    
    """
}

process abstract_gene_pool {
    publishDir "${params.outdir}/02.pre_gene/filtered_gene_info/", pattern: "*"
    tag "abstract_gene_pool"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(pool_vcf) from filtered_pool_ch
    output:
        tuple val(sample_name), path("${sample_name}.filtered_pool_gene.txt") into snpeff_pool_ch
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        ${snpsift} extractFields ${pool_vcf} "CHROM" "POS" "REF" "ALT" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" -s "\\t" -e "." > ${sample_name}.filtered_pool_gene.txt
    """
}

process anno_by_gnomAD {
    publishDir "${params.outdir}/02.pre_gene/anno_by_gnomAD/", pattern: "*"
    tag "anno_by_gnomAD"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(pool_vcf), path(tbi) from filtered_anno_vcf_ch1
    output:
        tuple val(sample_name), path("${sample_name}.filtered_pool_anno_by_gnomAD.vcf.gz") into anno_by_gnomAD_ch
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        ${snpsift} annotate ${gnomAD} ${pool_vcf} -info AC_eas,AF_eas,non_cancer_AC_eas,non_cancer_AF_eas > ${sample_name}.filtered_pool_anno_by_gnomAD.vcf.gz
    """
}

process anno_by_cosmic {
    publishDir "${params.outdir}/02.pre_gene/anno_by_cosmic/", pattern: "*"
    tag "anno_by_cosmic"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(pool_vcf), path(tbi) from filtered_anno_vcf_ch2
    output:
        tuple val(sample_name), path("${sample_name}.filtered_pool_anno_by_cosmic.vcf.gz") into anno_by_cosmic_ch
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        ${snpsift} annotate ${cosmic} ${pool_vcf} > ${sample_name}.filtered_pool_anno_by_cosmic.vcf.gz
    """
}

process anno_by_dbsnp {
    publishDir "${params.outdir}/02.pre_gene/anno_by_dbsnp/", pattern: "*"
    tag "anno_by_dbsnp"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(pool_vcf), path(tbi) from filtered_anno_vcf_ch3
        path change_chr_name_dbsnp from change_chr_name_dbsnp
    output:
        tuple val(sample_name), path("${sample_name}.filtered_pool_anno_by_dbsnp.vcf.gz") into anno_by_dbsnp_ch
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        cat ${change_chr_name_dbsnp} | awk '{print \$2,\$1}' > change_chr_name_dbsnp_reverse.txt
        bcftools annotate --rename-chrs ${change_chr_name_dbsnp} ${pool_vcf} -Oz -o new_pool_vcf.gz
        ${snpsift} annotate ${dbsnp} new_pool_vcf.gz | bcftools annotate --rename-chrs change_chr_name_dbsnp_reverse.txt -Oz -o ${sample_name}.filtered_pool_anno_by_dbsnp.vcf.gz
    """
}

merge_anno_ch = anno_by_gnomAD_ch
    .join(anno_by_dbsnp_ch, by : 0)
    .join(anno_by_cosmic_ch, by : 0)
    .join(snpeff_pool_ch, by : 0)

process filter_merge_anno {
    publishDir "${params.outdir}/02.pre_gene/filter_merge_anno/", pattern: "*"
    tag "filter_merge_anno"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(gnomAD_anno), path(dbsnp_anno), path(cosmic_anno), path(snpeff_anno) from merge_anno_ch
    output:
        tuple val(sample_name), path("${sample_name}.merge_anno.txt"), path("${sample_name}.merge_anno_for_dbNSFP.vcf") into merge_vcf_ch
        tuple val(sample_name), path("${sample_name}.merge_anno.txt"), path("${sample_name}.merge_anno_for_dbNSFP.vcf") into merge_vcf_ch2
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        bcftools filter -e "INFO/AF_eas > 0.0014" ${gnomAD_anno} -Oz -o ${sample_name}.gnomAD_filtered.vcf.gz
        bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF_eas\\n" ${sample_name}.gnomAD_filtered.vcf.gz > ${sample_name}.gnomAD_filtered.txt
        bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/COMMON\\n" ${dbsnp_anno} > ${sample_name}.dbsnp.txt
        bcftools filter -e "ID = '.'" ${cosmic_anno} -Oz -o ${sample_name}.cosmic_filtered.vcf.gz
        bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/HGVSP\\t%INFO/HGVSC\\n" ${sample_name}.cosmic_filtered.vcf.gz > ${sample_name}.cosmic_filtered.txt
        Rscript ${baseDir}/bin/merge_all_anno.R --cosmicfile ${sample_name}.cosmic_filtered.txt --dbsnpfile ${sample_name}.dbsnp.txt --gnomADfile ${sample_name}.gnomAD_filtered.txt --snpefffile ${snpeff_anno} --outmergefile ${sample_name}.merge_anno.txt --outvcffile ${sample_name}.merge_anno.vcf
        less ${sample_name}.gnomAD_filtered.vcf.gz | grep ^## > vcf.header.txt
        cat vcf.header.txt ${sample_name}.merge_anno.vcf > ${sample_name}.merge_anno_for_dbNSFP.vcf
    """
}
process anno_by_dbNSFP_hg38{
    publishDir "${params.outdir}/02.pre_gene/CPGs_result/", pattern: "*"
    tag "filter_merge_anno"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    input:
        tuple val(sample_name), path(merge_anno_txt), path(merge_anno_vcf) from merge_vcf_ch
    output:
        tuple val(sample_name), path("${sample_name}.CPGs_result.txt") into dbNSFP_result_ch
    script:
    """
        source  ~/app/bioinfo/dna/new.rc
        bcftools sort ${merge_anno_vcf} -o sort.vcf
        ${snpsift} dbnsfp -v -db ${dbNSFP_hg38} sort.vcf > anno_by_NFSP.vcf
        bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/dbNSFP_MetaSVM_pred\\t%INFO/dbNSFP_FATHMM_pred\\t%INFO/dbNSFP_LRT_pred\\t%INFO/dbNSFP_PROVEAN_pred\\t%INFO/dbNSFP_MutationTaster_pred\\t%INFO/dbNSFP_MutationAssessor_pred\\t%INFO/dbNSFP_Polyphen2_HDIV_pred\\t%INFO/dbNSFP_Polyphen2_HVAR_pred\\t%INFO/dbNSFP_SIFT_pred\\n" anno_by_NFSP.vcf > ${sample_name}.dbNSFP.txt
        Rscript ${baseDir}/bin/merge_dbSFSP.R --dbnfspfile ${sample_name}.dbNSFP.txt --mergeannofile ${merge_anno_txt} --outmergefile ${sample_name}.CPGs_result.txt
    """
}

//易感基因筛查
process pre_data{
	publishDir "${params.outdir}/02.pre_gene/CGC_anno" , pattern: "*"
	tag "pregene"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample_name), path(cpg) from dbNSFP_result_ch
    output:
        file "${sample_name}_CGC_anno.txt"
         file "${sample_name}_CGC_anno.xls"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    python3 ${baseDir}/bin/02.pregene_match.py ${cpg} /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hrda/ref/oncogene_database/cancer_gene_census.txt ${sample_name}_CGC_anno.txt 
    python3 ${baseDir}/bin/convert_txt_xls.py ${sample_name}_CGC_anno.txt 
    """    
}

//体细胞vcf解压
process somatic_vcf{
    publishDir "${params.outdir}/00.somatic_vcf/" , pattern: "*"
	tag "oncogene"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control) from pair_info3
    output:
        tuple val(sample),val(control),path("${sample}_${control}.vcf") into somatic_vcf
        tuple val(sample),val(control),path("${sample}_${control}.vcf") into somatic_vcf2
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    gunzip -c ${params.wgs_result_dir}/04.snpIndel_somatic/finalvcf/${sample}_vs_${control}_somatic_final_anno.vcf.gz > ${sample}_${control}.vcf
    """    
}
//驱动基因
process onco_data{
	publishDir "${params.outdir}/03.oncogene/" , pattern: "*"
	tag "oncogene2"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control),path(vcf) from somatic_vcf
    output:
        file "${sample}_${control}_oncogene.txt"
        file "${sample}_${control}_oncogene.xls"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    python3 ${baseDir}/bin/03.oncogene_match.py ${vcf} /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hrda/ref/oncogene_database/meiji_oncogene_database.txt ${sample}_${control}_oncogene.txt 
    python3 ${baseDir}/bin/convert_txt_xls.py ${sample}_${control}_oncogene.txt 
    """    
}
//药物靶点
process drug_target{
	publishDir "${params.outdir}/10.drug_target/result" , pattern: "*"
	tag "drugtarget"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control),path(vcf) from somatic_vcf2
    output:
        file "${sample}_${control}_drug_target.txt" into report_drug_ch
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    python3 ${baseDir}/bin/10.drug_result.py ${sample}_${control}.vcf /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hrda/ref/oncogene_database/meiji_drug_target.txt  temp.txt ${sample}_${control}
    awk 'NR==1{print;next} !a[\$2\$3]++' temp.txt | sort -k2,2V -k3,3n >> ${sample}_${control}_drug_target.txt
    python3 ${baseDir}/bin/convert_txt_xls.py ${sample}_${control}_drug_target.txt
    """    
}



//新抗原预测
process generate_somatic_vcf{
	publishDir "${params.outdir}//09.netMHC/full_somatic_vcf/" , pattern: "*"
	tag "generate_somatic_vcf"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "20G"
    input:
        tuple val(sample),val(control) from pair_info6
    output:
        tuple val(sample),val(control),path("${sample}_${control}_somatic_final_anno.vcf") into mhc_vcf
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    gunzip -c ${params.wgs_result_dir}/04.snpIndel_somatic/finalvcf/${sample}_vs_${control}_somatic_final_anno.vcf.gz > ${sample}_${control}.vcf
    python3 ${baseDir}/bin/09.vcf_filter.py ${sample}_${control}.vcf ${params.wgs_result_dir}/02.reference/ref.gtf ${sample}_${control}_somatic_final_anno.vcf
    """    
}

process generate_protein_fasta{
    publishDir "${params.outdir}/09.netMHC/protein_fasta" , pattern: "*"
    tag "mhc"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(sample),val(control),path(sample_all) from mhc_vcf
    output:
        tuple val(sample),val(control),file ("${sample}_${control}.protein.fasta") into protein_fasta,protein_fasta2
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    cp ${ref_protein} ref.protein.fa
    samtools faidx ref.protein.fa
    python3  ${baseDir}/bin/MHC_extract_protein.py ${sample_all} ref.protein.fa ./${sample}_${control}.protein.fasta
    """
}

mhc_I_ch = protein_fasta.combine(mhc_I).view()
mhc_II_ch = protein_fasta2.combine(mhc_II).view()

process antigenPrediction{
    publishDir "${params.outdir}/09.netMHC/MHC_I/" , pattern: "*"
    tag "mhc_I"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(sample),val(control),file(protein_fasta),val(anti) from mhc_I_ch
    output:
        tuple val(sample),val(control),val(anti),file("*.xls") into mhc_I_xls
        tuple val(sample),val(control),val(anti),file("*.txt") into mhc_I_txt
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/netMHCpan-4.1/netMHCpan -a ${anti} -f ${protein_fasta} -inptype 0 -xls -xlsfile ${sample}_${control}_${anti}_NetMHCpan_out.xls
        python3  ${baseDir}/bin/parse_outputs.py -f ${sample}_${control}_${anti}_NetMHCpan_out.xls  -w 2.0 -s 0.2 -o ${sample}_${control}_${anti}_MHC_out_parsed.txt
        sed -i 1d ${sample}_${control}_${anti}_MHC_out_parsed.txt
        """
}
process generate_gene_list{
    publishDir "${params.outdir}/09.netMHC/gene_list/" , pattern: "*"
    tag "mhc_I_result"
    queue "SANGERDEV"
    executor "slurm"
    cpus 4
    memory "10G"
    input:
        file(gtf) from gtf_ch
    output:
        file("gtf_gene.list") into gene_list,gene_list2
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        python3 ${baseDir}/bin/gtf_to_gene.py ${gtf} ./gtf_gene.list
        """
}

mhc_I_xls_gene = mhc_I_xls.combine(gene_list).view()

process antigenresult{
    publishDir "${params.outdir}/09.netMHC/MHC_I/final_xlsx" , pattern: "*"
    tag "mhc_I_xls"
    queue "SANGERDEV"
    executor "slurm"
     cpus 4
    memory "10G"
    input:
        tuple val(sample),val(control),val(anti),file(xls),file(gene_list) from mhc_I_xls_gene
    output:
        file("*")
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        python3 ${baseDir}/bin/map_mhc_xls_gene_new.py ${gene_list} ${xls} ./${sample}_${control}_${anti}_final
        """
}
process antigentxt{
    publishDir "${params.outdir}/09.netMHC/MHC_I/final_stat_txt" , pattern: "*"
    tag "mhc_I_txt"
    queue "SANGERDEV"
    executor "slurm"
     cpus 4
    memory "10G"
    input:
        tuple val(sample),val(control),val(anti),file(txt) from  mhc_I_txt.collect()
    output:
        file("*")
    script:
        def all_txt = txt.collect().join(" ")
        """
        source  ~/app/bioinfo/dna/new.rc
        cat ${baseDir}/bin/mhc_stat_header.txt ${all_txt} > ${sample}_${control}_mhc_I_final_stat.txt
        """
}
process antigenPrediction2{
    publishDir "${params.outdir}/09.netMHC/MHC_II/" , pattern: "*"
    tag "mhc_II"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(sample),val(control),file(protein_fasta),val(anti) from mhc_II_ch
    output:
        tuple val(sample),val(control),val(anti),file("*.xls") into mhc_II_xls
        tuple val(sample),val(control),val(anti),file("*.txt") into mhc_II_txt
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/netMHCIIpan-4.1/netMHCIIpan -a ${anti} -f ${protein_fasta} -inptype 0 -length 12 -xls -xlsfile ${sample}_${control}_${anti}_NetMHCIIpan_out.xls
        python3  ${baseDir}/bin/parse_outputs.py -f ${sample}_${control}_${anti}_NetMHCIIpan_out.xls  -w 2.0 -s 0.2 -o ${sample}_${control}_${anti}_MHCII_out_parsed.txt
        sed -i 1d ${sample}_${control}_${anti}_MHCII_out_parsed.txt
        """
}


mhc_II_xls_gene = mhc_II_xls.combine(gene_list2).view()

process antigenresult2{
    publishDir "${params.outdir}/09.netMHC/MHC_II/final_xlsx" , pattern: "*"
    tag "mhc_II_xls"
    queue "SANGERDEV"
    executor "slurm"
     cpus 4
    memory "10G"
    input:
        tuple val(sample),val(control),val(anti),file(xls),file(gene_list) from mhc_II_xls_gene
    output:
        file("*")
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        python3 ${baseDir}/bin/map_mhc_xls_gene_new.py ${gene_list} ${xls} ./${sample}_${control}_${anti}_final
        """
}
process antigentxt2{
    publishDir "${params.outdir}/09.netMHC/MHC_II/final_stat_txt" , pattern: "*"
    tag "mhc_II_txt"
    queue "SANGERDEV"
    executor "slurm"
     cpus 4
    memory "10G"
    input:
        tuple val(sample),val(control),val(anti),file(txt) from  mhc_II_txt.collect()
    output:
        file("*")
    script:
        def all_txt = txt.collect().join(" ")
        """
        source  ~/app/bioinfo/dna/new.rc
        cat ${baseDir}/bin/mhc_stat_header.txt ${all_txt} > ${sample}_${control}_mhc_II_final_stat.txt
        """
}
//高频突变

process convert_cram2bam_together{
    publishDir "${params.outdir}/07.significantly_mutated/bam/" , pattern: "*"
    tag "cram2bamall"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "150G"
    input:
        tuple val(sample),val(control) from pair_info8
    output:
        tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"),val(control), path("${control}.bam"), path("${control}.bam.bai")into bam_ch_all,bam_ch_all2,bam_ch_all3,bam_ch_all4
    script:
    """
    samtools view -bS ${params.wgs_result_dir}/01.fastq_qc/cram/${sample}.mkdup.cram -T ${ref} -@ 16 > ${sample}.origin.bam
    samtools sort ${sample}.origin.bam -@ 16 -o ${sample}.bam
    samtools index -b ${sample}.bam
    samtools view -bS ${params.wgs_result_dir}/01.fastq_qc/cram/${control}.mkdup.cram -T ${ref} -@ 16 > ${control}.origin.bam
    samtools sort ${control}.origin.bam -@ 16 -o ${control}.bam
    samtools index -b ${control}.bam
    """
}

process generate_bam_list{
    publishDir "${params.outdir}/07.significantly_mutated/bam_list/" , pattern: "*"
	tag "generate_bam_list"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "10G"
    input:
        tuple val(sample),val(control) from bam_ch_all.collect()
    output:
        file("bam.list") into bam_list
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    python3 ${baseDir}/bin/music2_bam_list.py ${params.sample_tumor_pair} ${params.outdir}/07.significantly_mutated/bam/
    echo finish
    """ 
}

process generate_maf_sig{
    publishDir "${params.outdir}/07.significantly_mutated/sig_maf/" , pattern: "*"
    tag "generatemaf"
    queue "SANGERDEV"
    executor "slurm"
    cpus 4
    memory "8G"
    input:
        tuple val(sample),val(control),path(full_vcf) from full_vcf3
    output:
        tuple val("${sample}_${control}"),path("${sample}_${control}.maf")into sig_maf_ch
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    bcftools view -s ${sample},${control} ${full_vcf} > ${sample}_${control}.vcf
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\tGT:DP:AD[\\t%GT:%DP:%AD]\\n' ${sample}_${control}.vcf -o ${sample}.${control}.vcf
    perl ${baseDir}/bin/snptomaf_for_sig.pl ${sample}.${control}.vcf 100 0.05
    """
}

process generate_roi{
    publishDir "${params.outdir}/07.significantly_mutated/roi/" , pattern: "*"
    tag "generateroi"
    queue "SANGERDEV"
    executor "slurm"
    cpus 2
    memory "4G"
    input:
        tuple val(name),file(mafs) from sig_maf_ch.collect()
    output:
        tuple(path("all.roi"),path("all.maf")) into roi_ch1, roi_ch2
        file("all.roi") into roi_ch
    script:
    def all_maf = mafs.collect().join(" ")
    """
    source  ~/app/bioinfo/dna/new.rc
    grep -v "sca" ${params.wgs_result_dir}/02.reference/ref.gtf |grep "exon"| cut -f1,4,5 > temp.1 && grep -v "sca"  ${params.wgs_result_dir}/02.reference/ref.gtf |grep "exon" | cut -d "\\"" -f12 > temp.2 && paste temp.1 temp.2 > all.roi
    cat ${baseDir}/bin/maf_header.txt ${all_maf}> all.maf
    
    """
}
process run_calc_covg_helper{
    publishDir "${params.outdir}/07.significantly_mutated/music2/roi_covgs/" , pattern: "*"
    tag "run_calc_covg_helper"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "16G"
    input:
        tuple val(roi), val(maf) from roi_ch1
        tuple val(sample), path(sample_bam), path(sample_bai),val(control), path(control_bam), path(control_bai) from bam_ch_all2
    output:
        path("${sample}.covg") into covg_ch
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        music2 bmr calc-covg-helper --normal-tumor-bam-pair \"${sample}:${control_bam}:${sample_bam}\" --roi-file ${roi} --reference-sequence ${ref} --output-file ./${sample}.covg --normal-min-depth 6 --tumor-min-depth 8 --min-mapq 20 --bp-class-types AT,CG,CpG
        """
}
process run_calc_covg{
    publishDir "${params.outdir}/07.significantly_mutated/music2/" , pattern: "*"
    tag "all"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "100G"
    input:
        tuple val(roi), val(maf) from roi_ch2
        file(covgs) from covg_ch.collect()
        file(bamlist) from bam_list
    output:
        file "*"
        file("smgs") into smg_ch
    script:
    def all_covg = covgs.collect().join(" ")
        """
        source  ~/app/bioinfo/dna/new.rc
        mkdir gene_covgs
        mkdir roi_covgs
        mkdir result
        Rscript ${baseDir}/bin/panorama_of_genomic_mutations.R --maf ${maf} --out ./result/
        cp -d ${all_covg} ./roi_covgs/
        sed "s/\t/:/g" ${bamlist} > bam_po.list
        music2 bmr calc-covg --roi-file ${roi} --reference-sequence ${ref} --bam-list  bam_po.list  --output-dir . 
        music2 bmr calc-bmr --roi-file ${roi} --reference-sequence ${ref} --bam-list ${bamlist} --maf-file ${maf} --output-dir . --show-skipped
        music2 smg --gene-mr-file gene_mrs --output-file smgs --max-fdr 0.05 --processors 1
        """
}

process orgdb{
    publishDir "${params.outdir}/07.significantly_mutated/enrich/", pattern:"*"
    executor 'slurm'
    queue "SANGERDEV"
    cpus 2
    cache 'lenient'
    input:
        file(smgs) from smg_ch
        file(roi) from roi_ch
    output:
        file "*"
    script:
        """
        source  ~/app/bioinfo/dna/new.rc
        mkdir LIB
        python3 ${baseDir}/bin/bsa_abstract_anno.py --infile ${anno_summary} --outfile go_kegg.list --go_anno GO_anno.csv --kegg_anno KEGG_rename.csv
        python3 ${baseDir}/bin/bsa_abstract_all_transcripts_from_anno_summary.py --infile ${anno_summary} --outfile genes_all.list
        Rscript ${baseDir}/bin/make_new_orgdb.R --annofile go_kegg.list --term2gene TERM2GENE.txt --term2name TERM2NAME.txt --outpath .
        R CMD INSTALL org.Ddemo.eg.db --library=./LIB/
        cut -f1 ${smgs} | sed '1d' > gene.list
        grep -f gene.list ${roi} > region.txt
        Rscript ${baseDir}/bin/abstract_pop_genes.R --regionfile region.txt --genefile genes_all.list --outfile genes_abstract.list
        cat genes_abstract.list | awk 'NR>1{print \$5}' > degfile.txt
        Rscript ${baseDir}/bin/enrich.R --degfile degfile.txt --term2genefile TERM2GENE.txt --term2namefile TERM2NAME.txt --outname all --db ./LIB/ --outdir all/ --title Significantly_high_frequency_Enrich_result
        """
}


//超突变样本
if (params.tumoronly == "no"){
process hla_bam{
    publishDir "${params.outdir}/06.hypermutation_classification/result/" , pattern: "*"
    tag "cram2bam"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "50G"
    input:
        tuple val(sample), path(sample_bam), path(sample_bai),val(control), path(control_bam), path(control_bai) from bam_ch_all3
    output:
        file "${sample}_${control}.prefix" into msi_report_ch
        file "${sample}_${control}.prefix_dis"
        file "${sample}_${control}.prefix_germline"
        file "${sample}_${control}.prefix_somatic"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    msisensor-pro scan -d ${ref} -o microsatellites.list
    msisensor-pro msi -d microsatellites.list -n ${control_bam} -t ${sample_bam} -o ${sample}_${control}.prefix
    """
}}else{
process msi_pro_tumor_only{
    publishDir "${params.outdir}/06.hypermutation_classification/tumor_only_result/" , pattern: "*"
    tag "msiprotumoronly"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "50G"
    input:
        tuple val(sample),val(control) from pair_info7
    output:
        file "${sample}_${control}.prefix" into msi_report_ch
        file "${sample}_${control}.prefix_all"
        file "${sample}_${control}.prefix_dis"
        file "${sample}_${control}.prefix_unstable"
    script:
    """
    source  ~/app/bioinfo/dna/new.rc
    samtools view -bS ${params.wgs_result_dir}/01.fastq_qc/cram/${sample}.mkdup.cram -T ${ref} -@ 8 > ${sample}.bam
    samtools sort ${sample}.bam -@ 8 -o ${sample}.sort.bam
    samtools index -b ${sample}.sort.bam
    msisensor-pro scan -d ${ref} -o microsatellites.list
    find ./ -name "${sample}.sort.bam" -printf "%p\t" > configure.txt
    msisensor-pro baseline -d microsatellites.list -i configure.txt -o ./
    msisensor-pro pro -d microsatellites.list -t ${sample}.sort.bam -o ${sample}_${control}.prefix
    """
}}
process msi_result{
        publishDir "${params.outdir}/06.hypermutation_classification/" , pattern: "*"
        tag "msi_result"
        queue "SANGERDEV"
        executor "slurm"
        cpus 4
        memory "10G"
        input:
            file (report) from msi_report_ch.collect() 
        output:
            file "*.xls"
            file "*.txt"
        script:
        def all_report = report.collect().join(" ")
            """
            source  ~/app/bioinfo/dna/new.rc
            echo ${all_report}
            python3 ${baseDir}/bin/05.hyper_result.py ./ ./ msi_pro
            """
}
//hla分型
process hla_subtype{
        publishDir "${params.outdir}/08.hla_subtype/sample/" , pattern: "*"
        tag "hla_subtype"
        queue "SANGERDEV"
        executor "slurm"
        cpus 8
        memory "50G"
        input:
            tuple val(sample), path(sample_bam), path(sample_bai),val(control), path(control_bam), path(control_bai) from bam_ch_all4
        output:
            file "*.report" into hla_report_ch
        script:
            """
            source  ~/app/bioinfo/dna/new.rc
            bash +e  ${baseDir}/bin/06.hla_scan.sh -b ${sample_bam} -s ${sample} -l ${baseDir}/bin/hla_list.txt
            """
}
process hla_result{
        publishDir "${params.outdir}/08.hla_subtype/result/" , pattern: "*"
        tag "hla_result"
        queue "SANGERDEV"
        executor "slurm"
        cpus 4
        memory "10G"
        input:
            file (report) from hla_report_ch.collect() 
        output:
            file "*.xls"
            file "*.txt"
        script:
        def all_report = report.collect().join(" ")
            """
            source  ~/app/bioinfo/dna/new.rc
            echo ${all_report}
            python3 ${baseDir}/bin/06.hla_report.py ./ ./
            """
}