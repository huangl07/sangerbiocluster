#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */
params.outdir = 'wgs_result'
params.help = false
params.bwameme = false
params.snpeff = false
params.cnvkit = false
params.gtf = 'None'
params.queue = 'SANGERDEV'
params.anno = 'None'
params.env = '~/app/bioinfo/dna/new.rc'
params.vcf = 'None'
params.target_bed = 'None'
params.vep_gtf = 'None'
params.script_dir = "${baseDir}"

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run genome_config.nf --outdir 'wgs_result'

    --ref        <file>  参考序列文件 [required]
    --anno       <file>  基因功能注释anno.detail
    --gtf        <file>  基因位置注释gtf
    --outdir     <dir>   publish文件夹
    --aslevel    <str>   参考基因组序列重命名文件 [required]
    --bwameme    <bool>  是否构建bwa-meme索引文件
    --snpeff     <bool>  是否构建snpeff索引
    --vcf        <file>  vcf文件，填了会自动运行snpeff
    --target_bed <file>  target_bed文件
    --cnvkit     <bool>  是否构建cnvkit索引
    --env        <file>  运行脚本前source的环境变量文件，默认为~/app/bioinfo/dna/new.rc
    --help       <bool>  打开帮助
    """.stripIndent()
}

def showParams() {
    log.info"""
    G E N O M E - C O N F I G  v 4.0
    ================================
    ref        : $params.ref
    anno       : $params.anno
    gtf        : $params.gtf
    outdir     : $params.outdir
    aslevel    : $params.aslevel
    bwameme    : $params.bwameme
    snpeff     : $params.snpeff
    vcf        : $params.vcf
    target_bed : $params.target_bed
    cnvkit     : $params.cnvkit
    env        : $params.env
    script_dir : $params.script_dir
    """.stripIndent()
}



process process_vepgtf {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 2
    memory '8G'
    queue "${params.queue}"
  input:
    path ref, name: 'raw.fa'
    path vepgtf, name: 'raw.gtf'
    path chromosome_list, name: 'raw_chromosome_list'

  output:
    path 'ref.vep.gtf.gz*', emit: gtf

  script:
    """
    source ${params.env}
    samtools faidx raw.fa
    python3 ${params.script_dir}/genome_config/check_chromosome_list.py -i raw_chromosome_list -f raw.fa.fai -o chromosome_list
    perl ${params.script_dir}/genome_config/GRename.pl -i raw.fa -o ref -f chromosome_list -g raw.gtf
    sort -k1,1V -k4,4n -k5,5n ref.gtf > sorted.gtf
    bgzip -c sorted.gtf > ref.vep.gtf.gz
    tabix -p gff ref.vep.gtf.gz
    """
}

process genome_grename {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 2
    memory '30G'
    queue "${params.queue}"
  input:
    path ref, name: 'raw.fa'
    val gtf
    path chromosome_list, name: 'raw_chromosome_list'

  output:
    path 'ref.fa', emit: ref
    path 'ref.gtf', emit: gtf
    path 'chromosome_list', emit: chlist
    path 'ref.genome.summary.xls', emit: summary
    path 'ref.table.xls', emit: table
    path 'ref.smtdict', emit: dict
    path 'ref.fa.fai', emit: fai
    path 'ref.stats.xls', emit: stats

  script:
    def gtfpar = gtf != 'None' ? '-g ' + gtf : ''
    """
  source ${params.env}
  samtools faidx raw.fa
  python3 ${params.script_dir}/genome_config/check_chromosome_list.py -i raw_chromosome_list -f raw.fa.fai -o chromosome_list
  perl ${params.script_dir}/genome_config/GRename.pl -i raw.fa -o ref -f chromosome_list ${gtfpar}
  seqkit fx2tab -j 4 -g -G -H -q -l -n ref.fa > ref.table.xls
  seqkit stats -a -T -j 4 ref.fa > ref.stats.xls
  python3 ${params.script_dir}/genome_config/merge_chrlist_stat.py -t ref.table.xls -c chromosome_list -o ref.genome.summary.xls
  samtools dict -o ref.smtdict ref.fa
  samtools faidx ref.fa
  """
}

process gff_check {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 2
    memory '20G'
    queue "${params.queue}"
  input:
    path gtf

  output:
    path 'gene_pos.txt', emit: pos

  script:
    """
    source ${params.env}
    python3 ${params.script_dir}/genome_config/gff_check.py --gtf ${gtf} --output .
    """
}

process rna_anno {
    publishDir "${params.outdir}/tmp/02.reference", pattern: 'anno.summary', mode: 'link'
    cpus 2
    memory '10G'
    queue "${params.queue}"
  input:
    path anno_detail
    path gene_pos

  output:
    path 'anno.summary', emit: anno

  script:
    """
    source ${params.env}
    python ${params.script_dir}/genome_config/get_anno_summary.py -i ${anno_detail} -o anno.summary -g ${gene_pos}
    """
}

process make_enrich_db {
    publishDir "${params.outdir}/tmp/02.reference/orgDB", mode: 'link', pattern: 'org.Ddemo.eg.db_1.0.tar.gz'
    publishDir "${params.outdir}/tmp/02.reference/orgDB", mode: 'link', pattern: 'kegg.gson'
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link', pattern: 'ref.transcripts.xls'
    cpus 2
    memory '24G'
    queue "${params.queue}"
  input:
    path anno

  output:
    path 'kegg.gson', emit: kegg
    path 'org.Ddemo.eg.db_1.0.tar.gz', emit: go
    path 'ref.transcripts.xls', emit: trans

  script:
    """
    source ${params.env}
    bash -e ${params.script_dir}/genome_config/anno_transcript.sh -a ${anno} -o ref.transcripts.xls
    Rscript ${params.script_dir}/genome_config/make_org_db.R -a ref.transcripts.xls -o .
    R CMD build org.Ddemo.eg.db
    """
}

process bgzip_tabix {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 4
    memory '8G'
    queue "${params.queue}"
  input:
    path ref
    val preset
    val out

  output:
    path '*.gz*', emit: out

  script:
    """
    source ${params.env}
    bgzip -c -@ 4 ${ref} > ${out}.gz
    tabix ${preset} ${out}.gz
    """
}

process cnvkit_index_wes {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 8
    memory '32G'
    queue "${params.queue}"
  input:
    path ref_fa, name: 'ref.fa'
    path ref_gtf, name: 'ref.gtf'
    path ref_bed, name: 'ref.bed'

  output:
    path 'ref.antitarget.bed'
    path 'ref.target.bed'
    path 'reference.cnn'

  script:
    """
    source ${params.env}
    cnvkit.py batch --process 8 --annotate ref.gtf --output-dir . -f ref.fa -n -t ref.bed
    """
}

process cnvkit_index_wgs {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 8
    memory '32G'
    queue "${params.queue}"
  input:
    path ref_fa, name: 'ref.fa'
    path ref_gtf, name: 'ref.gtf'

  output:
    path 'ref.antitarget.bed'
    path 'ref.target.bed'
    path 'reference.cnn'

  script:
    """
    source ${params.env}
    cnvkit.py batch --process 8 --annotate ref.gtf --output-dir . -f ref.fa -n -m wgs
    """
}

process bwa_index {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 16
    memory '128G'
    queue "${params.queue}"
  input:
    path ref , name: 'ref.fa'

  output:
    path 'ref.*', emit: out

  script:
    """
    source ${params.env}
    \$prefix/../miniconda3/envs/bwa-meme-2023/bin/bwa-meme index -a mem2 ref.fa -t 16
    \$prefix/../miniconda3/envs/bwa-meme-2023/bin/bwa-meme index -a meme ref.fa -t 16
    \$prefix/../miniconda3/envs/bwa-meme-2023/bin/bwa-meme-train-prmi -t 16 \
     --data-path \$(dirname ref.fa.suffixarray_uint64) ref.fa.suffixarray_uint64 \
     \$(basename ref.fa.suffixarray_uint64) pwl28,linear,linear_spline 268435456
    """
}

process gffread_to_fagtf {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 2
    memory '8G'
    queue "${params.queue}"
  input:
    path gff_path
    path fa

  output:
    path 'ref.cds.fa', emit: cds
    path 'ref.protein.fa', emit: prot
    path 'ref.mRNA.fa', emit: mrna

  script:
    """
    source ${params.env}
    gffread ${gff_path} -g ${fa} -x ref.cds.fa -w ref.mRNA.fa -y ref.protein.fa
    """
}

process snpeff_index {
    publishDir "${params.outdir}/tmp/02.reference", mode: 'link'
    cpus 2
    memory '60G'
    queue "${params.queue}"
  input:
    path ref_fa
    path ref_gff
    path ref_cds
    path ref_pro

  output:
    path 'snpEff.config', emit: conf
    path 'ref/*', emit: refdir

  script:
    """
    source ${params.env}
    rm -fr ref
    mkdir ref
    cp -l ${ref_fa} ref/sequences.fa
    cp -l ${ref_gff} ref/genes.gtf
    cp -l ${ref_cds} ref/cds.fa
    cp -l ${ref_pro} ref/protein.fa
    echo 'data.dir = .' > snpEff.config
    echo 'ref.genome : ref' >> snpEff.config
    java -Xmx10G -jar \$prefix/../snpEff/snpEff.jar build -gtf22 -v ref -c snpEff.config -noCheckCds -noCheckProtein
    """
}

process snpeff {
    publishDir "${params.outdir}/tmp/04.snpIndel/finalvcf", mode: 'link'
    cpus 8
    memory '60G'
    queue "${params.queue}"
  input:
    val ref_config
    path vcf, name: 'input.vcf'
    path 'chr.list'

  output:
    path 'pop.snpindel.final.vcf.gz*'

  script:
    """
    source ${params.env}
    awk -v OFS=" " '{print \$1,\$3}' chr.list > chr_rename.txt
    bcftools annotate --threads 8 --rename-chrs chr_rename.txt -Oz -o renamed.vcf.gz input.vcf
    java -Xmx10G -jar \$prefix/../snpEff/snpEff.jar eff -v ref -htmlStats snp.summary.html \
     -csvStats snpEff.snp.csv -c ${ref_config} renamed.vcf.gz > pop.snpindel.final.vcf
    bgzip -@ 8 pop.snpindel.final.vcf
    tabix pop.snpindel.final.vcf.gz || tabix -c pop.snpindel.final.vcf.gz
    """
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    } else {
        showParams()
    }
    if (params.vep_gtf != 'None'){
        process_vepgtf(file(params.ref), file(params.vep_gtf), file(params.aslevel))
    }
    if (params.gtf != 'None'){
        genome_grename(file(params.ref), file(params.gtf), file(params.aslevel))
        if (params.anno != 'None') {
            gff_check(genome_grename.out.gtf)
            rna_anno(file(params.anno), gff_check.out.pos)
            make_enrich_db(rna_anno.out.anno)
        }
        if (params.snpeff) {
            gffread_to_fagtf(genome_grename.out.gtf, genome_grename.out.ref)
            snpeff_index(
                genome_grename.out.ref,
                genome_grename.out.gtf,
                gffread_to_fagtf.out.cds,
                gffread_to_fagtf.out.prot
            )
            if (params.vcf != 'None') {
                snpeff(snpeff_index.out.conf, file(params.vcf), genome_grename.out.chlist)
            }
        }
        if (params.cnvkit) {
            if (params.target_bed != 'None'){
                cnvkit_index_wes(genome_grename.out.ref, genome_grename.out.gtf, file(params.target_bed))
        }else {
                cnvkit_index_wgs(genome_grename.out.ref, genome_grename.out.gtf)
            }
        }
    }else {
        genome_grename(file(params.ref), params.gtf, file(params.aslevel))
    }
    if (params.bwameme) {
        bwa_index(genome_grename.out.ref)
    }
    if (params.target_bed != 'None') {
        bgzip_tabix(file(params.target_bed), '-p bed', 'ref.bed')
    }
    workflow.onComplete {
        println 'G E N O M E - C O N F I G  v 4.0 Pipeline completed!'
    }
}


