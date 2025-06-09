#!/usr/bin/env nextflow
params.outdir = "primer_design_result"
params.help = false
params.env="~/app/bioinfo/dna/new.rc_1"
params.noregion=false
params.sites = 5000


def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    --vcf   <file>  input vcf file
    --outdir   <dir>   output dir
    --chr   <file>  chr list for draw
    ##############################################################
    --gff   <file>  gff file
    --anno_summary_file    <file>  anno_summary_file
    ##############################################################
    --group <file>   group file
        mbid must be given

        wpid    flag    phdep   pldep
        mpid    flag    phdep   pldep

    --popt  <str>   population
    ###############################################################
    --bulksize  <num>   bulk size
    --winsize <num>     windows size
    --stepsize <num>    step size
    ###############################################################
    --bootstrap <num>   bootstrap number
    --pvalue    <num>   pvalue
        ridit for pvalue [ 0.0001]
    ################################################################
    --grade <bool>  do grade for multi bulk than 2
    --mutmap <bool> do mutmap
    ###############################################################
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

if(!params.noregion){
    region=file(params.region)
    process extract_region_from_vcf{
        publishDir "${params.outdir}/01.extract_vcf", pattern: "extract.vcf.gz"
        queue "${params.queue}"
        executor "slurm"
        cpus 4
        memory "32G"
        input:
            file vcf_file from vcf
            file region_file from region
        output:
            file "extract.vcf" into extract_vcf_ch, extract_vcf_ch2
            file "extract.vcf.gz" into extract_gz_vcf_ch
        script:
            """
            source ${params.env}
            tabix -f ${vcf_file}
            bcftools view --threads 4 -R ${region_file} ${vcf_file} -o extract.vcf.gz
            gunzip -k extract.vcf.gz
            """
}
}else{
    process rename_gunzip_vcf{
        publishDir "${params.outdir}/01.extract_vcf", pattern: "extract.vcf.gz"
        queue "${params.queue}"
        executor "slurm"
        cpus 4
        memory "32G"
        input:
            file vcf_file from vcf
        output:
            file "extract.vcf" into extract_vcf_ch, extract_vcf_ch2
            file "extract.vcf.gz" into extract_gz_vcf_ch
        script:
            """
            source ${params.env}
            ln ${vcf_file} extract.vcf.gz
            file=\$(readlink -f extract.vcf.gz)
            rm extract.vcf.gz
            cp \$file extract.vcf.gz -rf
            tabix -f extract.vcf.gz
            gunzip -k extract.vcf.gz
            """
}
}

process primer3{
    publishDir "${params.outdir}/02.primer3", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    input:
        file extract_vcf_file from extract_vcf_ch
        file ref_file from ref
    output:
        file "variation.p3in"
        file "variation.p3out"
        file "variation.result"
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/bin:\$PATH
        export LD_LIBRARY_PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/lib64:\$LD_LIBRARY_PATH
        workdir=\$PWD
        perl ${baseDir}/bin/1.p3in.pl -d ${extract_vcf_file} -r 600-800 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./
        cd /mnt/lustre/users/sanger-dev/app/bioinfo/WGS/primer3/src
        ./primer3_core --output \$workdir/variation.p3out \$workdir/variation.p3in
        cd \$workdir
        perl ${baseDir}/bin/3.merge_p3out.v2.pl variation.p3out 3
        """
}

process split_vcf{
    publishDir "${params.outdir}/01.extract_vcf/split_vcf", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    input:
        path extract_vcf_file from extract_vcf_ch2
    output: 
        path "*" into split_vcf_ch
    script:
        """
        source ${params.env}
        java -jar /mnt/lustre/users/sanger-dev/app/bioinfo/dna/snpEff/SnpSift.jar split -l ${params.sites} ${extract_vcf_file}
        """
}

process filter_and_norm{
    publishDir "${params.outdir}/03.file_treatment/overlap_filter_norm", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    input:
        each split_vcf_file from split_vcf_ch
        file ref_file from ref
    output:
        tuple val("${key}"), path("${key}.overlap_filter_norm.vcf.gz") into norm_vcf_ch
    script:
        key=split_vcf_file.name.split("\\.")[1]
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        bgzip -c ${split_vcf_file} > split.vcf.gz
        tabix split.vcf.gz
        python3 ${baseDir}/bin/filter_overlap_pos.py --infile split.vcf.gz --outfile overlap_filter.list --no_overlap_len 100
        bcftools view -R overlap_filter.list split.vcf.gz -o overlap_filter.vcf.gz
        bcftools norm -f ${ref_file} -o ${key}.overlap_filter_norm.vcf.gz  -m - overlap_filter.vcf.gz
        """
}

process get_flank_seq{
    publishDir "${params.outdir}/03.file_treatment/t_polymarker_input", patetrn: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    input:
        tuple val(key), path(norm_vcf) from norm_vcf_ch
        file ref_file from ref
    output:
        file "${key}.flank_seq.fa" 
        tuple val(key), path("${key}.for_blast.fa"), path("${key}.t_polymarker_input.csv") into t_polymarker_input_ch
    script:
        """
        source ${params.env}
        python3 ${baseDir}/bin/get_pos.py ${norm_vcf} overlap_filter_norm.list
        python3 ${baseDir}/bin/get_flank_bed.py --infile overlap_filter_norm.list --outfile flank.bed --flank_len 50
        seqtk subseq ${ref_file} flank.bed > ${key}.flank_seq.fa
        python3 ${baseDir}/bin/get_input_config.py --infile ${key}.flank_seq.fa --outfile ${key}.t_polymarker_input.csv --snp_list overlap_filter_norm.list
        python3 ${baseDir}/bin/parse_polymarker_input.py ${key}.t_polymarker_input.csv ${key}.for_blast.fa
        """
}

process blast{
    publishDir "${params.outdir}/04.blast/blast_out", pattern: "*.blast_out.txt"
    publishDir "${params.outdir}/04.blast/flank_ranges", pattern: "*.flank_ranges.txt"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "32G"
    input:
        file ref_file from ref
        tuple val(key), path(for_blast_fa), path(t_polymarker_input) from t_polymarker_input_ch
    output:
        file "${key}.blast_out.txt"
        tuple val(key), path("${key}.flank_ranges.txt") into flank_ranges_ch, flank_ranges_ch2
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        makeblastdb -in ${ref_file} -dbtype nucl -parse_seqids
        blastn -task blastn -db ${ref_file} -query ${for_blast_fa} -outfmt "6 std qseq sseq slen" -word_size 11 -num_threads 8 -out ${key}.blast_out.txt
        python3 ${baseDir}/bin/getflanking.py ${t_polymarker_input} ${key}.blast_out.txt ${key}.flank_ranges.txt
        """
}

process kasp{
    publishDir "${params.outdir}/05.kasp/temp_marker", pattern: "temp_marker_*"
    publishDir "${params.outdir}/05.kasp", pattern: "KASP_output/selected_KASP_primers*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory "8G"
    input:
        file ref_file from ref
        tuple val(key), path(flank_ranges) from flank_ranges_ch
    output:
        file "temp_marker_*"
        file "KASP_output/selected_KASP_primers*" into kasp_result_ch
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        python3 ${baseDir}/bin/get_flank_bed2.py --infile ${flank_ranges}
        gawk '{ print \$2,\$3,\$4 > "temp_marker_"\$1".txt" }' ${flank_ranges}
        for i in flank_bed/*
        do
        base_name=`basename \${i}`
        seqtk subseq ${ref_file} \${i} > flanking_temp_marker_\${base_name}.fa
        done
        python3 /mnt/lustre/users/sanger-dev/app/bioinfo/dna/SNP_Primer_Pipeline2/bin/getkasp3.py 63 25 0
        """
}

process caps{
    publishDir "${params.outdir}/06.caps/temp_marker", pattern: "temp_marker_*"
    publishDir "${params.outdir}/06.caps", pattern: "CAPS_output/selected_CAPS_primers*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory "8G"
    input:
        file ref_file from ref
        tuple val(key), path(flank_ranges) from flank_ranges_ch2
    output:
        file "temp_marker_*"
        file "CAPS_output/selected_CAPS_primers*" into caps_result_ch
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        python3 ${baseDir}/bin/get_flank_bed2.py --infile ${flank_ranges}
        gawk '{ print \$2,\$3,\$4 > "temp_marker_"\$1".txt" }' ${flank_ranges}
        for i in flank_bed/*
        do
        base_name=`basename \${i}`
        seqtk subseq ${ref_file} \${i} > flanking_temp_marker_\${base_name}.fa
        done
        python3 /mnt/lustre/users/sanger-dev/app/bioinfo/dna/SNP_Primer_Pipeline2/bin/getCAPS.py 200 63 25 0
        """
}


process arrange_kasp_result{
    publishDir "${params.outdir}/05.kasp", pattern: "Potential_KASP_primers.tsv"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    input:
        file kasp_output from kasp_result_ch.collect()
    output:
        file "Potential_KASP_primers.tsv" into kasp_primers_ch
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        cat ${kasp_output} > Potential_KASP_primers.tsv
        """
}

process arrange_casp_result{
    publishDir "${params.outdir}/06.caps", pattern: "Potential_CAPS_primers.tsv"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    input:
        file caps_output from caps_result_ch.collect()
    output:
        file "Potential_CAPS_primers.tsv" into caps_primers_ch
    script:
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        cat ${caps_output} > Potential_CAPS_primers.tsv
        """
}

process result{
    publishDir "${params.outdir}/07.result", pattern: "caps_primers_result.xls"
    publishDir "${params.outdir}/07.result", pattern: "kasp_primers_result.xls"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    input:
        file kasp_primers from kasp_primers_ch
        file caps_primers from caps_primers_ch
    output:
        file "caps_primers_result.xls"
        file "kasp_primers_result.xls"
        """
        source ${params.env}
        export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH
        python3 ${baseDir}/bin/capskasp_trimlist.py ${caps_primers} caps_primers_result.xls ${kasp_primers} kasp_primers_result.xls
        """
}

