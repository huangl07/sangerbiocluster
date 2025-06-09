#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

params.norename = false
params.wgs = false
params.cog = false
params.cogped = ""
params.disease = ""
params.help = false
params.report = false
params.gene = ""
params.vcf = "published/data/04.snpIndel/finalvcf/pop.snpindel.final.vcf.gz"
params.noasso = false

def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run human_disease.nf --outdir hd_result --snpindel pop.vcf.gz --ped sample.ped

    --db        <file>  注释数据库所在文件夹，默认不写
    --wgsresult <file>  wgs output dir [required]
    --vcf       <file>  vcf相对wgs_result的相对路径，自定义vcf时才修改该参数
    --ped       <file>  样品信息表 [required]
    --disease   <file>  目标疾病名，一行一个
    --outdir    <dir>   输出文件夹
    --norename  <bool>  不需要重命名染色体
    --wgs       <bool>  是否wgs测序，默认wes
    --cog       <bool>  是否进行样品遗传模式分析
    --cogped    <file>  新生突变使用的家系文件
    --gene      <file>  老师感兴趣的基因列表
    --report    <bool>  是否生成报告
    --noasso    <bool>  不运行基因疾病关联分析
    --help      <bool>  打开帮助
    """.stripIndent()
}

def showParams(){
    log.info"""
    H U M A N - D I S E A S E  v 1.0
    ================================
    db        : $params.db
    wgsresult : $params.wgsresult
    vcf       : $params.vcf
    outdir    : $params.outdir
    ped       : $params.ped
    disease   : $params.disease
    norename  : $params.norename
    wgs       : $params.wgs
    cog       : $params.cog
    cogped    : $params.cogped
    gene      : $params.gene
    report    : $params.report
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
} else {
    showParams()
}

def is_set(option){
    if(option == ""){
        return ""
    } else {
        return file(option)
    }
}

outdir = file(params.outdir)
outdir.deleteDir()
result = outdir.mkdir()
println result ? "create directory: $outdir" : "Cannot create directory: $outdir"
params.outdir = file(params.outdir)

include { annot_snpindel } from "./module/snpindel"

process make_report {
    publishDir "${outdir}/publish", pattern: "*.pdf"
    publishDir "${outdir}/publish", pattern: "*.html"
    publishDir "${outdir}/publish", pattern: "data"
    cpus 2
    memory "100G"
    input:
        val finished
    output:
        path "data"
        path "*.pdf"
        path "*.html"
    script:
    def ltype = params.wgs ? "wgs" : "wes"
    def cog = params.cog ? 1 : 0
    def asso = params.noasso ? 0 : 1
    """
    cp -f ${projectDir}/bin/readme/结果说明文档.txt ${outdir}/annotation/
    bash -ue ${projectDir}/run_disease_report.sh ${params.wgsresult} ${outdir} . ${ltype} ${params.reportdb} ${cog} ${asso}
    """
}

workflow {
    if (params.wgs){
        type = "genome"
    } else {
        type = "exome"
    }
    disease = is_set(params.disease)
    cogped = is_set(params.cogped)
    gene = is_set(params.gene)
    annot_snpindel(
        params.wgsresult+'/'+ params.vcf,
        !params.norename,
        type,
        params.ped,
        params.cog,
        disease,
        cogped,
        gene,
        params.noasso
    )
    if (params.report){
        make_report(annot_snpindel.out.finished)
    }
}

workflow.onComplete {
    if(workflow.success){
        println "Workflow complete! ${workflow.duration}."
        println "HUMAN-DISEASE v 1.0 complete!"
    } else {
        println "Workflow failed! ${workflow.duration}."
    }
}