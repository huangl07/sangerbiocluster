#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

params.a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.kmer = 21
params.outdir = "gs_result"
params.help = false
params.report = false
params.name = "sample"
params.readlen = 151
params.project = ""
params.queue = "SANGERDEV"
params.basedir = "${baseDir}/bin"

def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run genome_survey.nf --outdir gs_result --name sample  --r1 r1.fq.gz --r2 r2.fq.gz --report -profile sangerdev

    --reportdb  <file>  报告模板文件夹 https://git.majorbio.com/ywtang/dna-report/-/tree/master/inst/rmarkdown/templates
    --r1        <file>  r1.fq.gz，required
    --r2        <file>  r2.fq.gz，required
    --outdir    <dir>   输出文件夹
    --name      <str>   样品名，required
    --a1        <str>   r1接头，外来数据设置为"Auto"
    --a2        <str>   r2接头，外来数据设置为"Auto"
    --kmer      <int>   kmer长度,奇数，默认21，基因组很大时可以酌情调大kmer
    --readlen   <int>   测序读长
    --report    <bool>  是否生成报告
    --project   <file>  project.info，required
    --help      <bool>  打开帮助
    """.stripIndent()
}

def showParams(){
    log.info"""
    G E N O M E - S U R V E Y v 1.0
    ================================
    reportdb  : $params.reportdb
    r1        : $params.r1
    r2        : $params.r2
    a1        : $params.a1
    a2        : $params.a2
    outdir    : $params.outdir
    name      : $params.name
    kmer      : $params.kmer
    readlen   : $params.readlen
    project   : $params.project
    report    : $params.report
    basedir   : $params.basedir
    """.stripIndent()
}





process fastp {
  publishDir "${params.outdir}", pattern: "*.fq.gz"
  publishDir "${params.outdir}/qc", pattern: "*.txt"
  publishDir "${params.outdir}/qc/json", pattern: "*.{xls,pdf,png}"
  cpus 8
  memory "32G"
  input:
    path r1
    path r2
    val name
    val a1
    val a2

  output:
    path "${name}_clean_1.fq.gz", emit: r1
    path "${name}_clean_2.fq.gz", emit: r2
    path "*.{txt,pdf,png,xls}", emit: results

  script:
  if(a1 == "Auto" || a2 == "Auto")
    """
    fastp -i ${r1} -I ${r2} -o ${name}_clean_1.fq.gz -O ${name}_clean_2.fq.gz -w 8 -h ${name}.html -j ${name}.json -M 20 -l 36 -n 0 -q 20 -3 20 -5 20
    mkdir -p qc
    perl ${params.basedir}/fastp.pl -i ${name}.json -o ${name}
    Rscript ${params.basedir}/ngsqc.r --base ${name}.raw.atgcn.xls --qual ${name}.raw.qual.xls --key ${name}.raw --od .
    Rscript ${params.basedir}/ngsqc.r --base ${name}.clean.atgcn.xls --qual ${name}.clean.qual.xls --key ${name}.clean --od .
    cat *.stat > qc.stat.txt
    """
  else
    """
    fastp -i ${r1} -I ${r2} -o ${name}_clean_1.fq.gz -O ${name}_clean_2.fq.gz -w 8 -h ${name}.html -j ${name}.json -M 20 -l 36 -n 0 -q 20 -3 20 -5 20  --adapter_sequence "${a1}" --adapter_sequence_r2 "${a2}"
    mkdir -p qc
    perl ${params.basedir}/fastp.pl -i ${name}.json -o ${name}
    Rscript ${params.basedir}/ngsqc.r --base ${name}.raw.atgcn.xls --qual ${name}.raw.qual.xls --key ${name}.raw --od .
    Rscript ${params.basedir}/ngsqc.r --base ${name}.clean.atgcn.xls --qual ${name}.clean.qual.xls --key ${name}.clean --od .
    cat *.stat > qc.stat.txt
    """
}

process smudgeplot {
    publishDir "${params.outdir}/smudgeplot", pattern: "*"
    cpus 8
    memory "64G"
    input:
        path r1
        path r2
        val name
        val kmer
    output:
        path "*.pdf", emit: pdf
        path "*.png", emit: png

    script:
    """
    mkdir -p FastK_Table
    FastK -v -t4 -T8 -k ${kmer} -M 64 -N FastK_Table ${r1} ${r2}
    smudgeplot.py hetmers -L 12 -t 8 -o kmerpairs --verbose FastK_Table
    smudgeplot.py plot -o ${name}_smu kmerpairs_text.smu
    """
}

process genome_survey {
  publishDir "${params.outdir}/gs", pattern: "*"
  cpus 8
  memory "64G"
  input:
    path r1
    path r2
    val name
    val kmer
    val readlen

  output:
    path "summary.txt", emit: summary
    path "linear_plot.png",emit: plot
    path "log_plot.png",emit: logplot
    val true, emit: finished

  script:
    """
    pigz -cd ${r1} > 1.fq
    pigz -cd ${r2} > 2.fq
    jellyfish-linux count -C -F 2 -m ${kmer} -s 1000000000 -t 8 -o ${name}.jf 1.fq 2.fq
    jellyfish-linux histo -t 8 ${name}.jf > ${name}.histo
    Rscript ${params.basedir}/genomescope2.R -i ${name}.histo -k ${kmer} -p 2 -o .
    """
}

process report  {
    publishDir "${params.outdir}/publish", pattern: "*.pdf"
    publishDir "${params.outdir}/publish", pattern: "*.html"
    publishDir "${params.outdir}/publish", pattern: "data"
    cpus 2
    memory "20G"
    input:
        val finished
        file outdir
        file project

    output:
        path "*.pdf", emit: pdf
        path "*.html", emit: html
        path "data", emit: data

    script:
    if(!params.report)
    """
    cp -f ${params.basedir}/readme/qc结果说明文档.txt ${outdir}/qc/结果说明文档.txt
    cp -f ${params.basedir}/readme/gs结果说明文档.txt ${params.outdir}/gs/结果说明文档.txt
    """
    else
    """
    cp -f ${params.basedir}/readme/qc结果说明文档.txt ${outdir}/qc/结果说明文档.txt
    cp -f ${params.basedir}/readme/gs结果说明文档.txt ${outdir}/gs/结果说明文档.txt
    bash -ue ${params.basedir}/run_report.sh ${outdir} . ${params.reportdb} ${project}
    """
}

workflow {
    if (params.help){
        helpMessage()
        exit 0
    } else {
        showParams()
    }
    fastp(file(params.r1),file(params.r2),params.name,params.a1,params.a2)
    genome_survey(fastp.out.r1, fastp.out.r2, params.name,params.kmer,params.readlen)
    report(genome_survey.out.finished, file(params.outdir), file(params.project))
}

workflow.onComplete {
    if(workflow.success){
        println "Workflow complete! ${workflow.duration}."
        println "GENOME-SURVEY v 1.0 complete!"
    } else {
        println "Workflow failed! ${workflow.duration}."
    }
}
