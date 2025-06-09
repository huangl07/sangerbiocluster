#!/usr/bin/env nextflow
nextflow.enable.dsl = 1
params.server="172"
params.outdir="insertion_check"
params.fq_list="fq.list"
params.queue="SANGERDEV"
params.env="~/app/bioinfo/dna/new.rc"
params.bwa=false
fq_list=Channel.from(file(params.fq_list))
              .splitCsv(header:false,sep:'\t')
params.ref="ref.fa"
params.insert="insert.txt"
adapter = "/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/adapter.fa"
ref_ch=Channel.from(file(params.ref))
ref_ch2=Channel.from(file(params.ref))
insert_ch=Channel.from(file(params.insert))
insert_ch1=Channel.from(file(params.insert))
if (params.server == "172" ){
    params.env="~/app/bioinfo/dna/new.rc"
    params.queue="SANGERDEV"
    conda_env = "~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh && set -h && conda activate bwa-meme-2023"
    sentieon = "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/sentieon-genomics-202112.07/bin/"
    bwa_meme = "~/app/bioinfo/dna/miniconda3/envs/bwa-meme-2023/bin/bwa-meme"
    bwa_meme2 = "~/app/bioinfo/dna/env/bin/bwa-mem2"
}else{
    params.env="/mnt/ilustre/users/xiaoya.ye/demo/BSA/1"
    conda_env = "~dna/.bash_conda"
    params.queue="DNA"
    sentieon = "/mnt/ilustre/users/dna/.env/sentieon-genomics-202112.06/bin/sentieon"
    bwa_meme = "/mnt/ilustre/users/dna/.env/BWA-MEME/bwa-meme"
    bwa_meme2 = "/mnt/ilustre/users/dna/.env/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
}

process fastp{
    publishDir "${params.outdir}/fq" , pattern: "*"
    tag "fastp"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "20G"
    input:
        tuple val(fq1),val(fq2),val(name) from fq_list
    output:
        tuple file("${name}_clean_1.fq.gz"),file("${name}_clean_2.fq.gz"),val(name) into clean_fq
        tuple file("${name}_clean_1.fq.gz"),file("${name}_clean_2.fq.gz"),val(name) into clean_fq1
        tuple file("${name}_clean_1.fq.gz"),file("${name}_clean_2.fq.gz"),val(name) into clean_fq2
        tuple file("${name}_clean_1.fq.gz"),file("${name}_clean_2.fq.gz"),val(name) into clean_fq3
        tuple file("${name}_clean_1.fq.gz"),file("${name}_clean_2.fq.gz"),val(name) into clean_fq4
        file("*")
    script:

    """
    source ${params.env}
    fastp -i ${fq1} -I ${fq2} -o ${name}_clean_1.fq.gz -O ${name}_clean_2.fq.gz -w 8 -h ${name}.html -j ${name}.json -M 20 -l 150 -n 0 -q 20 -3 20 -5 20
    perl ${baseDir}/bin/fastp.pl -i ${name}.json -o ${name}
    Rscript ${baseDir}/bin/ngsqc.r --base ${name}.raw.atgcn --qual ${name}.raw.qual --key ${name}.raw --od ./
    Rscript ${baseDir}/bin/ngsqc.r --base ${name}.clean.atgcn --qual ${name}.clean.qual --key ${name}.clean --od ./
    """
}

fq_list_ref = clean_fq4.combine(ref_ch2).view()
if (params.server == "172" ){
    process get_cram{
        publishDir "${params.outdir}/cram" , pattern: "*"
        tag "fastp"
        queue "${params.queue}"
        executor "slurm"
        cpus 8
        memory "64G"
        input:
            tuple val(fq1),val(fq2),val(sample),file(ref) from fq_list_ref
        output:
            tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch
            tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch1
            tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch2
        script:
        """
        source ${params.env}
        source ${conda_env}
        cp ${ref} ref.fa || echo ${ref}
        bwa-meme index -a mem2 ref.fa -t 8
        bwa-meme index -a meme ref.fa -t 8
        build_rmis_dna.sh ref.fa
        samtools faidx  ref.fa
        bwa-meme mem -7 -M -a -t 8 -R "@RG\\tID:${sample}\\tLG:${sample}\\tLB:1\\tPL:illumina\\tSM:${sample}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ref.fa ${fq1} ${fq2} > temp.sam
        samtools sort --reference ref.fa -o ${sample}.cram -O CRAM -@ 8 temp.sam
        rm -rf temp.sam
        samtools index -@ 8 ${sample}.cram
        mkdir ${sample}_temp
        export SENTIEON_LICENSE=/mnt/lustre/users/sanger-dev/app/bioinfo/WGS/sentieon-genomics-201808/MajorBio_cluster_0.37.lic
        export PATH=/mnt/lustre/users/sanger-dev/app/program/Python/bin:/bin:/usr/bin:/opt/rocks/bin:\$PATH
        export PYTHONPATH=/mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/src:/mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/webroot:\$PYTHONPATH
        export LANG=en_US.UTF-8
        sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo LocusCollector --fun score_info ${sample}.score.txt
        sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo Dedup --score_info ${sample}.score.txt  --metric ${sample}.metric.txt ${sample}.mkdup.cram
        """
    }
}else{
    process get_cram_236{
    publishDir "${params.outdir}/cram" , pattern: "*"
    tag "fastp"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "64G"
    input:
        tuple val(fq1),val(fq2),val(sample),file(ref) from fq_list_ref
    output:
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch1
        tuple val(sample),path("${sample}.mkdup.cram"),path("ref.fa"),path("${sample}.mkdup.cram.crai"),path("${sample}.mkdup.cram.bai") into cram_ch2
    script:
    """
    source ${params.env}
    source ${conda_env}
    cp ${ref} ref.fa || echo ${ref}
    bwa-meme index -a mem2 ref.fa -t 8
    bwa-meme index -a meme ref.fa -t 8
    build_rmis_dna.sh ref.fa
    samtools faidx  ref.fa
    bwa-meme mem -7 -M -a -t 8 -R "@RG\\tID:${sample}\\tLG:${sample}\\tLB:1\\tPL:illumina\\tSM:${sample}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ref.fa ${fq1} ${fq2} > temp.sam
    samtools sort --reference ref.fa -o ${sample}.cram -O CRAM -@ 8 temp.sam
    rm -rf temp.sam
    samtools index -@ 8 ${sample}.cram
    mkdir ${sample}_temp
    sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo LocusCollector --fun score_info ${sample}.score.txt
    sentieon driver --temp_dir ./${sample}_temp -t 8 -i ${sample}.cram --reference ref.fa --algo Dedup --score_info ${sample}.score.txt  --metric ${sample}.metric.txt ${sample}.mkdup.cram
    """
}
}

process samtoolstat{
    publishDir "${params.outdir}/samtools_stat" , pattern: "*"
    tag "fastp"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "64G"
    input:
        tuple val(sample),path(sample_cram),path(ref),path(crai),path(bai) from cram_ch
    output:
        file("*")
    script:
    """
    source ${params.env}
    samtools stats -@ 8 --reference ${ref} ${sample_cram} > ${sample}.mapstat
    perl ${baseDir}/bin/map_stat.pl -b ${sample}.mapstat -i ${sample}.insert.xls -c ${sample}.depth.xls -o ${sample}.summary.stat -k ${sample} -f ${ref}.fai
    Rscript ${baseDir}/bin/insert_size.R -i ${sample}.insert.xls -o ${sample}.insert -s ${sample}
    Rscript ${baseDir}/bin/coverage_depth.R -i ${sample}.depth.xls -o ${sample}.depth -d ${sample}
    """
}

process mosdepthstat{
    publishDir "${params.outdir}/mosdepth_stat" , pattern: "*"
    tag "fastp"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "64G"
    input:
        tuple val(sample),path(sample_cram),path(ref),path(crai),path(bai) from cram_ch1
    output:
        file("*")
    script:
    """
    source ${params.env}
    mosdepth -t 8 -f ${ref} --fast-mode -b 100000 ${sample} ${sample_cram}
    Rscript ${baseDir}/bin/genomeCoveragehorizontalArea.R --infile ${sample}.regions.bed.gz  --idfile ${ref}.fai --outfile ${sample}.genome.coverage --group.col 1 --x.col 2 --y.col 4 --x.lab Sequence-Position --y.lab AverageDepth-log2 --title.lab 'genome coverage on ${sample}' --skip 0 --unit kb --log2
    """
}

process tdna{
    publishDir "${params.outdir}/tdna/" , pattern: "*"
    tag "tdna"
    queue "${params.queue}"
    executor "slurm"
    cpus 32
    memory "500G"
    input:
        tuple val(fq1),val(fq2),val(name)from clean_fq
    output:
        file("*")
    script:
    """
    source ${params.env}
    python3 ${baseDir}/bin/tdnascan.py -1 ${fq1} -2 ${fq2} -g ${params.ref} -t ${params.insert} -p ${name} || echo 'Failed!!' > ${name}.tdna.status
    """
}

process aimhii{
    publishDir "${params.outdir}/aimhii/" , pattern: "*"
    tag "aimhii"
    queue "${params.queue}"
    executor "slurm"
    cpus 32
    memory "200G"
    input:
        tuple val(fq1),val(fq2),val(name) from clean_fq1
    output:
        file("*")
    script:
    """
    source ${params.env}
    mkdir ${name}_temp
    aimhii ${params.ref} ${params.insert} ${adapter} ${fq1} ${fq2} --outfile ${name}.csv --plot ${name} -t `readlink -f ${name}_temp` --threads 16 || echo 'Failed!!' > ${name}.aimhii.status
    """
}

process ref_set{
    publishDir "${params.outdir}/ref/" , pattern: "*"
    tag "refset"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "64G"
    input:
        path(ref) from ref_ch
        path(insert) from insert_ch
    output:
        tuple path("insert_ref.fa"),path("insert_name") into index
        //tuple path("ref.fa"),path("ref.fa.0123"),path("ref.fa.amb"),path("ref.fa.ann"),path("ref.fa.bwt.2bit.64"),path("ref.fa.pac") into index
        //tuple path("insert_ref.fa"),path("insert_ref.fa.0123"),path("insert_ref.fa.amb"),path("insert_ref.fa.ann"),path("insert_ref.fa.bwt.2bit.64"),path("insert_ref.fa.pac"),path("insert_name")into bwa
    script:
    if(params.bwa)
    """
    source ${params.env}
    cat ${ref} ${insert} > insert_ref.fa
    bwa index insert_ref.fa
    grep ">" ${insert} > insert_name
    sed -i "s/>//g" insert_name
    """
    else
    """
    source ${params.env}
    cat ${ref} ${insert} > insert_ref.fa
    bwa-mem2 index insert_ref.fa
    grep ">" ${insert} > insert_name
    sed -i "s/>//g" insert_name
    """
}

//ref_index_ch = clean_fq3.combine(index).view()
ref_bwa_ch = clean_fq2.combine(index).view()
//ins_ref_index_ch = ref_index_ch.combine(insert_ch1).view()

process ref_cov{
    publishDir "${params.outdir}/ref_cov/" , pattern: "*"
    tag "refcov"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "50G"
    input:
        tuple val(fq1),val(fq2),val(name),path(insert_fa),path(insert_name) from ref_bwa_ch
    output:
        file("${name}.depth")
        file("${name}*.png")
    script:
    if(params.bwa)
    """
    source ${params.env}
    bwa mem \$(readlink -f ${insert_fa}) ${fq1} ${fq2} -t 8 > ${name}.sam
    samtools view -b ${name}.sam | samtools sort -o ${name}.sorted.bam
    samtools index ${$name}.sorted.bam
    samtools faidx ${insert_fa}
    less -S insert_ref.fa.fai | grep -f insert_name | perl -ne '@a=split;print "samtools depth -r \$a[0]:1-\$a[1] ${name}.sorted.bam\n"' | bash > ${name}.depth
    Rscript ${baseDir}/bin/draw_depth.R --file ${name}.depth --name ${name} --insert_name ${insert_name}
    rm -rf ${name}.sam
    """
    else
    """
    source ${params.env}
    bwa-mem2 mem ${insert_fa} ${fq1} ${fq2} -t 8 > ${name}.sam
    samtools view -b ${name}.sam | samtools sort -o ${name}.sorted.bam
    samtools index ${$name}.sorted.bam
    samtools faidx ${insert_fa}
    less -S insert_ref.fa.fai | grep -f insert_name | perl -ne '@a=split;print "samtools depth -r \$a[0]:1-\$a[1] ${name}.sorted.bam\n"' | bash > ${name}.depth
    Rscript ${baseDir}/bin/draw_depth.R --file ${name}.depth --name ${name} --insert_name ${insert_name}
    rm -rf ${name}.sam
    """
}

//process TLOC {
//	publishDir "${params.outdir}/tloc", pattern: "*"
//	tag "TLOC"
//    queue "${params.queue}"
//    cpus 16
//    memory "64G"
//    executor "slurm"
//    input:
//        tuple val(fq1),val(fq2),val(name),path(index),path(index_123),path(index_amb),path(index_ann),path(index_bit),path(index_pac),path(insert) from ins_ref_index_ch
//    output:
//        file ("*")
//    script:
//    """
//    source ${params.env}
//    python2.7 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/demo/output/test_py/T-LOC-mem/T-LOC2.py  --fastq ${fq1},${fq2} --genome_Bwa ${index}  --TDNA ${insert} --output ./${name} --Sample_name ${name} --genome ${index}
//    """
//}
