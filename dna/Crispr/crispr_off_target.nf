#!/usr/bin/env nextflow
params.help = false
params.outdir = "demo"
params.env="~/app/bioinfo/dna/new.rc_1"
params.queue="SANGERDEV"
def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:


    --vcf   <file>  input vcf file
    --ref   <file>  input ref fasta file
    --outdir   <dir>   output dir
    --chr   <file>  chr list for draw
    ##############################################################
    --gff   <file>  gff file
    --pam   <file>   pam file
    --sgrna <file>  sgrna file
    ##############################################################
    --group <file>   group file
        sample2    case    100  5
        sample1    control    100   5
    ###############################################################
    --case_id       case sample id
    --control_id    control sample id
    ##############################################################
    --queue <str> slurm partion
        dna for 236
        SANGERDEV for 172
    --env   <str>   env for some package
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
vcf=file(params.vcf)
ref=file(params.ref)
chr=file(params.chr)
pam=file(params.pam)
sgrna=file(params.sgrna)
case_id=params.case_id
control_id=params.control_id
group_list=file(params.group_list)

process split_fa{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/00.pretreatment", pattern:"ref"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory '8GB'
    input:
        file ref from ref
        file chr from chr
    output:
        path "ref" into ref_ch
    script:
        """
        faidx -x ${ref}
        mkdir ref
        for i in `cat ${chr}`
        do
        mv \${i}.fa ref
        done
        """
}

process index_genome{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/00.pretreatment", pattern:"genome_library"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory '8G'
    input:
        path ref from ref_ch
        file pam from pam
    output:
        path "genome_library" into index_genome_ch
    script:
        """
        source /mnt/lustre/users/sanger-dev/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
        conda activate crispritz
        python3.8 ~/app/bioinfo/dna/miniconda3/envs/crispritz/bin/crispritz.py index-genome ref_fa ${ref} ${pam} -bMax 2 -th 4
        conda deactivate
        """
}

process search_homo_region{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/01.homo_region", pattern:"*.xls"
    publishDir "${params.outdir}/01.homo_region", pattern:"*.txt"
    publishDir "${params.outdir}/03.result", pattern:"*.pdf"
    publishDir "${params.outdir}/03.result", pattern:"*.png"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory '8G'
    input:
        path index_genome from index_genome_ch
        path ref from ref_ch
        file pam from pam
        file sgrna from sgrna
    output:
        path "homo_region.extended_profile.xls"
        path "homo_region.profile.xls"
        path "homo_region.profile_complete.xls"
        path "homo_region.profile_dna.xls"
        path "homo_region.profile_rna.xls"
        path "homo_region.scores.txt"
        path "homo_region.targets.CFD.txt" , optional: true
        path "homo_region.targets.txt" into homo_region_ch
        path "homo_region_seqlogo.pdf"
        path "homo_region_seqlogo.png"
    script:
        """
        source /mnt/lustre/users/sanger-dev/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
        conda activate crispritz
        python3.8 ~/app/bioinfo/dna/miniconda3/envs/crispritz/bin/crispritz.py search ${index_genome}/* ${pam} ${sgrna} homo_region -index -mm 5 -bDNA 2 -bRNA 2 -t -scores ${ref}
        Rscript ${baseDir}/bin/seqlogo.R --infile homo_region.extended_profile.xls --outfile1 homo_region_seqlogo.pdf --outfile2 homo_region_seqlogo.png
        conda deactivate
        """
}

process arrange_and_get_bed{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/01.homo_region", pattern:"homo_region.bed.xls"
    publishDir "${params.outdir}/03.result", pattern:"homo_region.result.xls"
    publishDir "${params.outdir}/03.result", pattern:"homo_region.stat.xls"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory '8G'
    input:
        path homo_region from homo_region_ch
    output:
        path "homo_region.result.xls"
        path "homo_region.bed.xls" into homo_region_bed_ch
        path "homo_region.stat.xls"
    script:
        """
        python3 ${baseDir}/bin/arrange_homo_region.py --infile ${homo_region} --outfile homo_region.result.xls --bedfile homo_region.bed.xls
        python3 ${baseDir}/bin/stat_homo_region.py --infile homo_region.result.xls --outfile homo_region.stat.xls
        """
}

process extract_homo_region_from_vcf{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/02.filter_vcf", pattern:"homo_region.extract.vcf.gz"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory '8G'
    input:
        file vcf from vcf
        file homo_region_bed from homo_region_bed_ch
        val case_id
        val control_id
    output:
        path "homo_region.extract.vcf.gz" into extract_vcf_ch
    script:
        """
        tabix -f ${vcf}
        #bcftools view --threads 4 -s ${case_id},${control_id} -R ${homo_region_bed} ${vcf} -o temp.vcf.gz
        #bcftools annotate -a ${homo_region_bed} -c
        bcftools view --threads 4 -s ${case_id},${control_id} -R ${homo_region_bed} ${vcf} -o homo_region.extract.vcf.gz
        """
}

process filter_vcf{
    beforeScript "source ${params.env}"
    publishDir "${params.outdir}/02.filter_vcf", pattern:"vcf_table.xls"
    publishDir "${params.outdir}/03.result", pattern:"targe_info.xls"
    publishDir "${params.outdir}/03.result", pattern:"targe_stat.xls"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory '8G'
    input:
        path vcf from extract_vcf_ch
        path group_list from group_list
    output:
        path "vcf_table.xls"
        path "targe_info.xls"
        path "targe_stat.xls"
    script:
        """
        python3 ${baseDir}/bin/vcf2table.py --infile ${vcf} --outfile vcf_table.xls
		grep -v "None" vcf_table.xls > vcf_table2.xls
        python3 ${baseDir}/bin/table_filter.py --table_file vcf_table2.xls --group_file ${group_list} --lowdepth 10 --outfile filter_table.xls
        python3 ${baseDir}/bin/split_table.py --vcf_table filter_table.xls --outfile targe_info.xls
        python3 ${baseDir}/bin/stat_table.py --infile targe_info.xls --outfile targe_stat.xls
        """
}


workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
