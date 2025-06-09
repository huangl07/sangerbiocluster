#!/usr/bin/env nextflow
params.sample_tumor_pair="sample_tumor_pair"
params.wgs_result_dir="wgs_result_dir"
params.hla="yes"
params.hlalist="HLA-A"
params.outdir="./"
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
process convert_cram2bam{
    publishDir "${params.outdir}/05.hypermutation_classification/bam_sample/" , pattern: "*"
    tag "cram2bam"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "16G"
    input:
        tuple val(sample),val(control) from pair_info
    output:
        tuple val(sample), file("${sample}.sort.bam"), file("${sample}.sort.bam.bai")into bam_ch1,bam_ch1_2
    script:
    """
    samtools view -bS ${params.wgs_result_dir}/01.fastq_qc/cram/${sample}.mkdup.cram -T ${params.wgs_result_dir}/genome_config_output/ref.fa -@ 8 > ${sample}.bam
    samtools sort ${sample}.bam -@ 8 -o ${sample}.sort.bam
    samtools index -b ${sample}.sort.bam
    """
}
process convert_cram2bam_control{
    publishDir "${params.outdir}/05.hypermutation_classification/bam_control/" , pattern: "*"
    tag "cram2bam"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "16G"
    input:
        tuple val(sample),val(control) from pair_info2
    output:
        tuple val(control), file("${control}.sort.bam"), file("${control}.sort.bam.bai")into bam_ch2
    script:
    """
    samtools view -bS ${params.wgs_result_dir}/01.fastq_qc/cram/${control}.mkdup.cram -T ${params.wgs_result_dir}/genome_config_output/ref.fa -@ 8 > ${control}.bam
    samtools sort ${control}.bam -@ 8 -o ${control}.sort.bam
    samtools index -b ${control}.sort.bam
    """
}
process hypermutation_classification{
        publishDir "${params.outdir}/05.hypermutation_classification/sample/" , pattern: "*"
        tag "hypermutation_classification"
        queue "SANGERDEV"
        executor "slurm"
        cpus 8
        memory "50G"
        input:
            tuple val(sample),sample_bam,sample_bam_bai from bam_ch1
        output:
            file "*"
        script:
            """
            source  ~/app/bioinfo/dna/new.rc
            /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/msi/msisensor2/msisensor2 msi -b 8 -M /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/msi/msisensor2/models_hg38/ -t ${sample_bam} -o ./${sample}.tumor.prefix
            """
}
process hypermutation_classification_control{
        publishDir "${params.outdir}/05.hypermutation_classification/control/" , pattern: "*"
        tag "hypermutation_classification"
        queue "SANGERDEV"
        executor "slurm"
        cpus 8
        memory "50G"
        input:
            tuple val(control),control_bam,control_bam_bai from bam_ch2
        output:
            file "*"
        script:
            """
            source  ~/app/bioinfo/dna/new.rc
            /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/msi/msisensor2/msisensor2 msi -b 8 -M /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/msi/msisensor2/models_hg38/ -t ${control_bam} -o ./${control}.tumor.prefix
            """
}
process hla_subtype{
        publishDir "${params.outdir}/07.hla_subtype/sample/" , pattern: "*"
        tag "hla_subtype"
        executor "local"
        input:
            tuple val(sample),sample_bam,sample_bam_bai from bam_ch1_2
        output:
            file  "${sample}.hlatype.report" into report
            file "*"
        script:
            """
            source  ~/app/bioinfo/dna/new.rc
            /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/hla_scan_r_v2.1.4 -v 38 -d /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/db/HLA-ALL.IMGT -g ${params.hlalist} -b ${sample_bam} -t 8 -f 75 > ${sample}.hlatype.report 2>${sample}.error
            if [ $? == 1 ] then
                touch ${sample}.noresult.log
            fi
            """
    }
/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/hla_scan_r_v2.1.4 -v 38 -d /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/db/HLA-ALL.IMGT -g ${params.hlalist} -b ${sample_bam} -t 8 -f 75 > ${sample}.hlatype.report  1>${sample}.noresult.log 2>${sample}.error