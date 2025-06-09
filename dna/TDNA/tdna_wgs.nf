#!/usr/bin/env nextflow

params.outdir="insertion_check"
params.fq_list="fq.list"
fq_list=Channel.from(file(params.fq_list))
              .splitCsv(header:false,sep:'\t')
fq_list1=Channel.from(file(params.fq_list))
              .splitCsv(header:false,sep:'\t')
fq_list2=Channel.from(file(params.fq_list))
              .splitCsv(header:false,sep:'\t')
fq_list3=Channel.from(file(params.fq_list))
              .splitCsv(header:false,sep:'\t')
params.ref="ref.fa"
params.insert="insert.txt"
adapter = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/adapter.fa"
ref_ch=Channel.from(file(params.ref))
insert_ch=Channel.from(file(params.insert))
insert_ch1=Channel.from(file(params.insert))

process tdna{
    publishDir "${params.outdir}/tdna/" , pattern: "*"
    tag "tdna"
    queue "SANGERDEV"
    executor "slurm"
    cpus 32
    memory "500G"
    input:
        tuple val(fq1),val(fq2),val(name)from fq_list
    output:
        file("*")
    script:
    """
    source  /mnt/lustre/users/sanger-dev/app/bioinfo/dna/new.rc
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
        tuple val(fq1),val(fq2),val(name) from fq_list1
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

ref_index_ch = fq_list2.combine(index).view()
ref_bwa_ch = fq_list3.combine(bwa).view()
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
//    memory "200G"
//    executor "slurm"
//    input:
//        tuple val(fq1),val(fq2),val(name),path(index),path(index_123),path(index_amb),path(index_ann),path(index_bit),path(index_pac),path(insert) from ins_ref_index_ch
//    output:
//        file("*")
//    script:
//    """
//    source ~/app/bioinfo/dna/new.rc
//    python2.7 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/demo/output/test_py/T-LOC-mem/T-LOC2.py  --fastq ${fq1},${fq2} --genome_Bwa ${index}  --TDNA ${insert} --output ./${name} --Sample_name ./${name} --genome ${index}
//    """
//}
