#!/usr/bin/env nextflow


/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */


params.outdir="tdna_result"
params.fq_list = "fq.list"
params.a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.help = false
params.project = "project.txt"
params.report = false
params.insert = "insert.fa"
params.ref = "ref.fa"

def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run tdna_v2.nf --report -profile sangerdev

    Arguments:
    -profile    <str>  The profile to use, 'sangerdev'(172 default) / 'dna'(236)
    --reportdb  <file>  报告模板文件夹 https://git.majorbio.com/ywtang/dna-report/-/tree/master/inst/rmarkdown/templates
    --fq_list   <file>  fq list，default ./fq.list
    --outdir    <dir>   输出文件夹，default ./tdna_result
    --a1        <str>   r1接头，华大测序时可以为""
    --a2        <str>   r2接头，华大测序时可以为""
    --insert    <file>  插入片段文件，default ./insert.fa
    --ref       <file>  参考基因组文件，default ./ref.fa
    --report    <bool>  是否生成报告， default no
    --project   <file>  project.info文件，default ./project.txt
    --help      <bool>  打开帮助
    """.stripIndent()
}



def showParams(){
    log.info"""
    G E N O M E - S U R V E Y v 1.0
    ================================
    reportdb  : $params.reportdb
    fq_list   : $params.fq_list
    a1        : $params.a1
    a2        : $params.a2
    outdir    : $params.outdir
    insert    : $params.insert
    ref       : $params.ref
    project   : $params.project
    report    : $params.report
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
} else {
    showParams()
}

outdir = file(params.outdir)
outdir.deleteDir()
result = outdir.mkdir()
println result ? "create directory: $outdir" : "Cannot create directory: $outdir"
outdir = file(params.outdir)

process fastp {
    publishDir "${params.outdir}/qc" , pattern: "*"
    tag "${name}"
    cpus 16
    memory "20G"
    input:
        tupe val(name),path(fq1),path(fq2)
        val a1
        val a2
    output:
        tuple val(name), path("${name}_clean_1.fq.gz"), path("${name}_clean_2.fq.gz"), emit: fq
        path "*.{txt,pdf,png,xls}", emit: results
        path "*.stat", emit: stat
    script:

    """
    source ~/app/bioinfo/dna/new.rc
    fastp -i ${r1} -I ${r2} -o ${name}_clean_1.fq.gz -O ${name}_clean_2.fq.gz -w 8 -h ${name}.html -j ${name}.json -M 20 -l 36 -n 0 -q 20 -3 20 -5 20  --adapter_sequence "${a1}" --adapter_sequence_r2 "${a2}"
    perl ${projectDir}/bin/fastp.pl -i ${name}.json -o ${name}
    Rscript ${projectDir}/bin/ngsqc.r --base ${name}.raw.atgcn.xls --qual ${name}.raw.qual.xls --key ${name}.raw --od .
    Rscript ${projectDir}/bin/ngsqc.r --base ${name}.clean.atgcn.xls --qual ${name}.clean.qual.xls --key ${name}.clean --od .
    """
}

process qc_stat {
    publishDir "${params.outdir}/qc" , pattern: "*"
    cpus 1
    memory "1G"
    input:
       path '*'
    output:
        path "qc.stat.txt"
    script:
        """
        cat *.stat > qc.stat.txt
        """
}

process ref_index {
    publishDir "${params.outdir}/reference" , pattern: "*"

}

process get_cram{
    publishDir "${params.outdir}/cram" , pattern: "*"
    tag "${name}"
    cpus 8
    memory "64G"
    input:
        tuple val(fq1),val(fq2),val(sample),file(ref)from fq_list_ref
    output:
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch1
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch2
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    source ~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
    set -h
    conda activate bwa-meme-2023
    cp ${ref} ref.fa || echo ${ref}
    ~/app/bioinfo/dna/miniconda3/envs/bwa-meme-2023/bin/bwa-meme index -a mem2 ref.fa -t 8
    ~/app/bioinfo/dna/miniconda3/envs/bwa-meme-2023/bin/bwa-meme index -a meme ref.fa -t 8
    build_rmis_dna.sh ref.fa
    ~/app/bioinfo/dna/env/bin/samtools faidx  ref.fa
    ~/app/bioinfo/dna/miniconda3/envs/bwa-meme-2023/bin/bwa-meme mem -7 -M -a -t 8 -R "@RG\\tID:${sample}\\tLG:${sample}\\tLB:1\\tPL:illumina\\tSM:${sample}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ref.fa ${fq1} ${fq2} > temp.sam
    ~/app/bioinfo/dna/env/bin/samtools sort --reference ref.fa -o ${sample}.cram -O CRAM -@ 8 temp.sam
    rm -rf temp.sam
    ~/app/bioinfo/dna/env/bin/samtools index -@ 8 ${sample}.cram
    mkdir ${sample}_temp
    export SENTIEON_LICENSE=/mnt/lustre/users/sanger-dev/app/bioinfo/WGS/sentieon-genomics-201808/MajorBio_cluster_0.37.lic
    export PATH=/mnt/lustre/users/sanger-dev/app/program/Python/bin:/bin:/usr/bin:/opt/rocks/bin:\$PATH
    export PYTHONPATH=/mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src:/mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/webroot:\$PYTHONPATH
    export LANG=en_US.UTF-8
    /mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/sentieon-genomics-202112.07/bin/sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo LocusCollector --fun score_info ${sample}.score.txt
    /mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/sentieon-genomics-202112.07/bin/sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo Dedup --score_info ${sample}.score.txt  --metric ${sample}.metric.txt ${sample}.mkdup.cram
    """
}

process samtoolstat{
    publishDir "${params.outdir}/samtools_stat" , pattern: "*"
    tag "fastp"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "64G"
    input:
        tuple val(sample),path(sample_cram),path(ref),path(crai),path(bai) from cram_ch
    output:
        file("*")
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    ~/app/bioinfo/dna/env/bin/samtools stats -@ 8 --reference ${ref} ${sample_cram} > ${sample}.mapstat
    ~/app/bioinfo/dna/env/bin/perl /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src/mbio/packages/wgs_v4/map_stat.pl -b ${sample}.mapstat -i ${sample}.insert.xls -c ${sample}.depth.xls -o ${sample}.summary.stat -k ${sample} -f ${ref}.fai
    ~/app/bioinfo/dna/env/bin/Rscript /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src/mbio/packages/wgs_v4/insert_size.R -i ${sample}.insert.xls -o ${sample}.insert -s ${sample}
    /mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/Rscript /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src/mbio/packages/wgs_v4/coverage_depth.R -i ${sample}.depth.xls -o ${sample}.depth -d ${sample}
    """
}

process mosdepthstat{
    publishDir "${params.outdir}/mosdepth_stat" , pattern: "*"
    tag "fastp"
    queue "SANGERDEV"
    executor "slurm"
    cpus 8
    memory "64G"
    input:
        tuple val(sample),path(sample_cram),path(ref),path(crai),path(bai) from cram_ch1
    output:
        file("*")
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    /mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/mosdepth -t 8 -f ${ref} --fast-mode -b 100000 ${sample} ${sample_cram}
    /mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/Rscript /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src/mbio/packages/wgs_v4/genomeCoveragehorizontalArea.R --infile ${sample}.regions.bed.gz  --idfile ${ref}.fai --outfile ${sample}.genome.coverage --group.col 1 --x.col 2 --y.col 4 --x.lab Sequence-Position --y.lab AverageDepth-log2 --title.lab 'genome coverage on ${sample}' --skip 0 --unit kb --log2
    """
}

process tdna{
    publishDir "${params.outdir}/tdna/" , pattern: "*"
    tag "tdna"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "500G"
    input:
        tuple val(fq1),val(fq2),val(name)from clean_fq
    output:
        file("*")
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    ${baseDir}/bin/tdnascan.py -1 ${fq1} -2 ${fq2} -g ${params.ref} -t ${params.insert} -p ${name}
    """
}

process aimhii{
    publishDir "${params.outdir}/aimhii/" , pattern: "*"
    tag "aimhii"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "200G"
    input:
        tuple val(fq1),val(fq2),val(name) from clean_fq1
    output:
        file("*")
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    mkdir ${name}_temp
    aimhii ${params.ref} ${params.insert} ${adapter} ${fq1} ${fq2} --outfile ${name}.csv --plot ${name} -t `readlink -f ${name}_temp` --threads 16
    """
}

process ref_set{
    publishDir "${params.outdir}/ref/" , pattern: "*"
    tag "refset"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "64G"
    input:
        path(ref) from ref_ch
        path(insert) from insert_ch
    output:
        path "ref.fa*"
        tuple path("ref.fa"),path("ref.fa.0123"),path("ref.fa.amb"),path("ref.fa.ann"),path("ref.fa.bwt.2bit.64"),path("ref.fa.pac") into index
        tuple path("insert_ref.fa"),path("insert_ref.fa.0123"),path("insert_ref.fa.amb"),path("insert_ref.fa.ann"),path("insert_ref.fa.bwt.2bit.64"),path("insert_ref.fa.pac"),path("insert_name")into bwa
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    cp ${ref} ref.fa || echo ${ref}
    ~/app/bioinfo/dna/env/bin/bwa-mem2 index ref.fa
    cat ${ref} ${insert} > insert_ref.fa
    ~/app/bioinfo/dna/env/bin/bwa-mem2 index insert_ref.fa
    grep ">" ${insert} > insert_name
    sed -i "s/>//g" insert_name
    """
}

ref_index_ch = clean_fq3.combine(index).view()
ref_bwa_ch = clean_fq2.combine(bwa).view()
ins_ref_index_ch = ref_index_ch.combine(insert_ch1).view()

process ref_cov{
    publishDir "${params.outdir}/ref_cov/" , pattern: "*"
    tag "refcov"
    queue "SANGERDEV"
    executor "slurm"
    cpus 16
    memory "50G"
    input:
        tuple val(fq1),val(fq2),val(name),path(insert_fa),path(insert_fa_123),path(insert_fa_amb),path(insert_fa_ann),path(insert_fa_bit),path(insert_fa_pac),path(insert_name) from ref_bwa_ch
    output:
        file("${name}.depth")
        file("${name}*.png")
    script:
    """
    source ~/app/bioinfo/dna/new.rc
    ~/app/bioinfo/dna/env/bin/bwa-mem2 mem ${insert_fa} ${fq1} ${fq2} -t 8 > ${name}.sam
    samtools view -b ${name}.sam | samtools sort -o ${name}.sorted.bam
    samtools depth ${name}.sorted.bam > ${name}.depth
    Rscript ${baseDir}/bin/draw_depth.R --file ${name}.depth --name ${name} --insert_name ${insert_name}
    """
}

//process TLOC {
//	publishDir "${params.outdir}/tloc", pattern: "*"
//	tag "TLOC"
//    queue "SANGERDEV"
//    cpus 16
//    memory "64G"
//    executor "slurm"
//    input:
//        tuple val(fq1),val(fq2),val(name),path(index),path(index_123),path(index_amb),path(index_ann),path(index_bit),path(index_pac),path(insert) from ins_ref_index_ch
//    output:
//        file ("*")
//    script:
//    """
//    source ~/app/bioinfo/dna/new.rc
//    python2.7 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/demo/output/test_py/T-LOC-mem/T-LOC2.py  --fastq ${fq1},${fq2} --genome_Bwa ${index}  --TDNA ${insert} --output ./${name} --Sample_name ${name} --genome ${index}
//    """
//}
