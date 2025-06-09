#!/usr/bin/env nextflow
params.outdir = "primer_design_result"
params.help = false
params.env="~/app/bioinfo/dna/new.rc_1"
params.queue="SANGERDEV"
params.snp_only=false
params.indel_only=false

def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow-1 run -bg primer_design_multi.nf --vcf {filtered_vcf} --ref {ref} --region_list {region}

    --vcf           <file>     vcf.gz [required]
    --raw           <file>     pop.snpindel.final.vcf.gz [required]
    --region_list   <file>     region.list文件 [required]
    --snp_only      <bool>     不对indel设计
    --env           <file>     source的环境 默认 "~/app/bioinfo/dna/new.rc_1"
    --queue         <string>   slurm队列 默认 "SANGERDEV"
    --outdir        <dir>      输出文件夹 默认 "primer_design_result"
    --help          <bool>     打开帮助
    """.stripIndent()
}

def showParams(){
    log.info"""
    G E N O M E - C O N F I G  v 4.0
    ================================
    vcf          : $params.vcf
    raw          : $params.raw
    ref          : $params.ref
    region_list  : $params.region_list
    snp_only     : $params.snp_only
    indel_only   : $params.indel_only
    queue        : $params.queue
    outdir       : $params.outdir
    env          : $params.env
    help         : $params.help
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
} else {
    showParams()
}

if (params.queue != "SANGERDEV" ){
    p3_path = "/mnt/ilustre/users/xiaoya.ye/BSA/primer3/src"
    env_path = "/mnt/ilustre/users/xiaoya.ye/BSA"
    export_path = ""
    source_env = ""
}else{
    p3_path = "/mnt/lustre/users/sanger-dev/app/bioinfo/WGS/primer3/src"
    env_path = "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/"
    export_path = "export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/align/ncbi-blast-2.12.0+/bin:\$PATH"
    source_env = "source ${params.env}"
}

vcf=file(params.vcf)
ref=file(params.ref)
regions = Channel.from(file(params.region_list))
              .splitCsv(header:false,sep:'\t')
              .groupTuple()
raw=file(params.raw)

process extract_region_from_vcf{
    publishDir "${params.outdir}/01.extract_vcf", pattern: "${method}.extract.vcf.gz"
    publishDir "${params.outdir}/01.extract_vcf", pattern: "${method}.extract.snp.vcf.gz"
    publishDir "${params.outdir}/01.extract_vcf", pattern: "${method}.extract.indel.vcf.gz"
    queue "${params.queue}"
    executor "slurm"
    tag "${method}"
    cpus 4
    memory "32G"
    input:
        file vcf_file from vcf
        tuple val(method), path(region_file) from regions
    output:
        tuple val(method), path("${method}.extract.vcf") into extract_vcf_ch, extract_vcf_ch2
        tuple val(method), path("${method}.extract.snp.vcf") into extract_snp_vcf_ch
        tuple val(method), path("${method}.extract.indel.vcf") into extract_indel_vcf_ch
        tuple "${method}.extract.vcf.gz"
    script:
        """
        ${source_env}
        tabix -f ${vcf_file} || tabix -f -C ${vcf_file}
        bcftools view --threads 4 -R ${region_file} ${vcf_file} -o ${method}.extract.vcf.gz
        bcftools view --threads 4 --types snps ${method}.extract.vcf.gz -o ${method}.extract.snp.vcf
        bcftools view --threads 4 --types indels ${method}.extract.vcf.gz -o ${method}.extract.indel.vcf
        gunzip -k ${method}.extract.vcf.gz
        """
}

if (!params.indel_only && params.queue == "SANGERDEV" ){
process primer3_snp{
    publishDir "${params.outdir}/02.primer3", pattern: "${method}_snp"
    publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.sanger.snp.result.xls"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "20G"
    tag "${method}_snp"
    input:
        tuple val(method), path(extract_vcf_file) from extract_snp_vcf_ch
        file ref_file from ref
    output:
        path "${method}_snp"
        path "${method}.sanger.snp.result.xls"
    script:
        """
        ${source_env}
        mkdir ${method}_snp
        export PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/bin:\$PATH
        export LD_LIBRARY_PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/lib64:\$LD_LIBRARY_PATH
        workdir=\$PWD
        perl ${baseDir}/bin/1.p3in.pl -d ${extract_vcf_file} -r 600-800 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./${method}_snp
        cd ${p3_path}
        ./primer3_core --output \$workdir/${method}_snp/variation.p3out \$workdir/${method}_snp/variation.p3in
        cd \$workdir
        perl ${baseDir}/bin/3.merge_p3out.v2.pl \$workdir/${method}_snp/variation.p3out 3
        python3 ${baseDir}/bin/primer_result_arrange.py --infile \$workdir/${method}_snp/variation.result --outfile ${method}.tmp.sanger.snp.result.xls
        less -S ${method}.tmp.sanger.snp.result.xls|grep -v Chrom|perl -ne '@a=split;\$a[12]=reverse(\$a[12]);\$a[12]=~tr/ATGC/TACG/;print ">",join(",",\$a[0],\$a[1],\$a[-2],\$a[-1]),"\\n",\$a[8].("X" x 10).\$a[12],"\\n"' > primer.fa
        samtools faidx ref.fa
         less -S ref.fa.fai |perl -ne '@a=split;print "samtools faidx ref.fa \$a[0]  > \$a[0].fa && blat  -tileSize=11 -stepSize=5 -t=DNA  -q=DNA \$a[0].fa primer.fa sub.\$a[0].fa.psl  -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead\\n";'|bash
         cat sub.*.psl > primer.fa.psl
        #blat -tileSize=11 -stepSize=5 -t=DNA  -q=DNA  ref.fa primer.fa primer.fa.psl -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead
        less -S primer.fa.psl |perl -ne '{@a=split;\$stat{\$a[9]}++ if(\$a[0]/(\$a[10]-10) > 0.95)}END{foreach \$a(keys %stat){print \$a if(\$stat{\$a} > 1)}}' > multi.list
        less -S ${method}.tmp.sanger.snp.result.xls multi.list|perl -ne 'if(/Chrom/){print \$_}else{chomp;@a=split;if(scalar @a > 4){\$id=join(",",\$a[0],\$a[1],\$a[-2],\$a[-1])}else{\$id=\$_};\$stat{\$id}++;\$info{\$id}=\$_;if(\$stat{\$id} == 1){print \$info{\$id},"\\n"}}' >  ${method}.sanger.snp.result.xls
        """
}
}else if (!params.indel_only && params.queue != "SANGERDEV") {
    process primer3_snp_236{
    publishDir "${params.outdir}/02.primer3", pattern: "${method}_snp"
    publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.sanger.snp.result.xls"
    queue "${params.queue}"
    executor "slurm"
    tag "${method}_snp"
    cpus 1
    memory "20G"
    input:
        tuple val(method), path(extract_vcf_file) from extract_snp_vcf_ch
        file ref_file from ref
    output:
        path "${method}_snp"
        path "${method}.sanger.snp.result.xls"
    script:
        """
        #${source_env}
        mkdir ${method}_snp
        workdir=\$PWD
        perl ${baseDir}/bin/1.p3in.pl -d ${extract_vcf_file} -r 600-800 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./${method}_snp
        cd ${p3_path}
        ./primer3_core --output \$workdir/${method}_snp/variation.p3out \$workdir/${method}_snp/variation.p3in
        cd \$workdir
        perl ${baseDir}/bin/3.merge_p3out.v2.pl \$workdir/${method}_snp/variation.p3out 3
        python3 ${baseDir}/bin/primer_result_arrange.py --infile \$workdir/${method}_snp/variation.result --outfile ${method}.tmp.sanger.snp.result.xls
        less -S ${method}.tmp.sanger.snp.result.xls|grep -v Chrom|perl -ne '@a=split;\$a[12]=reverse(\$a[12]);\$a[12]=~tr/ATGC/TACG/;print ">",join(",",\$a[0],\$a[1],\$a[-2],\$a[-1]),"\\n",\$a[8].("X" x 10).\$a[12],"\\n"' > primer.fa
            samtools faidx ref.fa
         less -S ref.fa.fai |perl -ne '@a=split;print "samtools faidx ref.fa \$a[0]  > \$a[0].fa && blat  -tileSize=11 -stepSize=5 -t=DNA  -q=DNA \$a[0].fa primer.fa sub.\$a[0].fa.psl  -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead\\n";'|bash
         cat sub.*.psl > primer.fa.psl
        #blat -tileSize=11 -stepSize=5 -t=DNA  -q=DNA  ref.fa primer.fa primer.fa.psl -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead
        less -S primer.fa.psl |perl -ne '{@a=split;\$stat{\$a[9]}++ if(\$a[0]/(\$a[10]-10) > 0.95)}END{foreach \$a(keys %stat){print \$a if(\$stat{\$a} > 1)}}' > multi.list
        less -S ${method}.tmp.sanger.snp.result.xls multi.list|perl -ne 'if(/Chrom/){print \$_}else{chomp;@a=split;if(scalar @a > 4){\$id=join(",",\$a[0],\$a[1],\$a[-2],\$a[-1])}else{\$id=\$_};\$stat{\$id}++;\$info{\$id}=\$_;if(\$stat{\$id} == 1){print \$info{\$id},"\\n"}}' >  ${method}.sanger.snp.result.xls
        """
}
}
if(!params.snp_only && params.queue == "SANGERDEV" ){
    process primer3_indel{
        publishDir "${params.outdir}/02.primer3", pattern: "${method}_indel"
        publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.sanger.indel.result.xls"
        queue "${params.queue}"
        executor "slurm"
        tag "${method}_indel"
        cpus 1
        memory "20G"
        input:
            tuple val(method), path(extract_vcf_file) from extract_indel_vcf_ch
            file ref_file from ref
        output:
            path "${method}_indel"
            path "${method}.sanger.indel.result.xls"
        script:
            """
            ${source_env}
            mkdir ${method}_indel
            export PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/bin:\$PATH
            export LD_LIBRARY_PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/lib64:\$LD_LIBRARY_PATH
            workdir=\$PWD
            perl ${baseDir}/bin/1.p3in.pl -d ${extract_vcf_file} -r 80-150 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./${method}_indel
            cd ${p3_path}
            ./primer3_core --output \$workdir/${method}_indel/variation.p3out \$workdir/${method}_indel/variation.p3in
            cd \$workdir
            perl ${baseDir}/bin/3.merge_p3out.v2.pl \$workdir/${method}_indel/variation.p3out 3
            python3 ${baseDir}/bin/primer_result_arrange.py --infile \$workdir/${method}_indel/variation.result --outfile ${method}.tmp.sanger.indel.result.xls
            less -S ${method}.tmp.sanger.indel.result.xls|grep -v Chrom|perl -ne '@a=split;\$a[12]=reverse(\$a[12]);\$a[12]=~tr/ATGC/TACG/;print ">",join(",",\$a[0],\$a[1],\$a[-2],\$a[-1]),"\\n",\$a[8].("X" x 10).\$a[12],"\\n"' > primer.fa
                samtools faidx ref.fa
         less -S ref.fa.fai |perl -ne '@a=split;print "samtools faidx ref.fa \$a[0]  > \$a[0].fa && blat  -tileSize=11 -stepSize=5 -t=DNA  -q=DNA \$a[0].fa primer.fa sub.\$a[0].fa.psl  -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead\\n";'|bash
         cat sub.*.psl > primer.fa.psl
            #blat -tileSize=11 -stepSize=5 -t=DNA  -q=DNA  ref.fa primer.fa primer.fa.psl -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead
            less -S primer.fa.psl |perl -ne '{@a=split;\$stat{\$a[9]}++ if(\$a[0]/(\$a[10]-10) > 0.95)}END{foreach \$a(keys %stat){print \$a if(\$stat{\$a}>1)}}' > multi.list
            less -S ${method}.tmp.sanger.indel.result.xls multi.list|perl -ne 'if(/Chrom/){print \$_}else{chomp;@a=split;if(scalar @a > 4){\$id=join(",",\$a[0],\$a[1],\$a[-2],\$a[-1])}else{\$id=\$_};\$stat{\$id}++;\$info{\$id}=\$_;if(\$stat{\$id} == 1){print \$info{\$id},"\\n"}}' >  ${method}.sanger.indel.result.xls
            """
    }
}else if ( !params.snp_only && params.queue != "SANGERDEV" ){
        process primer3_indel_236{
        publishDir "${params.outdir}/02.primer3", pattern: "${method}_indel"
        publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.sanger.indel.result.xls"
        queue "${params.queue}"
        executor "slurm"
        tag "${method}_indel"
        cpus 1
        memory "20G"
        input:
            tuple val(method), path(extract_vcf_file) from extract_indel_vcf_ch
            file ref_file from ref
        output:
            path "${method}_indel"
            path "${method}.sanger.indel.result.xls"
        script:
            """
            mkdir ${method}_indel
            workdir=\$PWD
            perl ${baseDir}/bin/1.p3in.pl -d ${extract_vcf_file} -r 80-150 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./${method}_indel
            cd ${p3_path}
            ./primer3_core --output \$workdir/${method}_indel/variation.p3out \$workdir/${method}_indel/variation.p3in
            cd \$workdir
            perl ${baseDir}/bin/3.merge_p3out.v2.pl \$workdir/${method}_indel/variation.p3out 3
            python3 ${baseDir}/bin/primer_result_arrange.py --infile \$workdir/${method}_indel/variation.result --outfile ${method}.tmp.sanger.indel.result.xls
            less -S ${method}.tmp.sanger.indel.result.xls|grep -v Chrom|perl -ne '@a=split;\$a[12]=reverse(\$a[12]);\$a[12]=~tr/ATGC/TACG/;print ">",join(",",\$a[0],\$a[1],\$a[-2],\$a[-1]),"\\n",\$a[8].("X" x 10).\$a[12],"\\n"' > primer.fa
               samtools faidx ref.fa
         less -S ref.fa.fai |perl -ne '@a=split;print "samtools faidx ref.fa \$a[0]  > \$a[0].fa && blat  -tileSize=11 -stepSize=5 -t=DNA  -q=DNA \$a[0].fa primer.fa sub.\$a[0].fa.psl  -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead\\n";'|bash
         cat sub.*.psl > primer.fa.psl
           # blat -tileSize=11 -stepSize=5 -t=DNA  -q=DNA  ref.fa primer.fa primer.fa.psl -minScore=10 -minIdentity=95 -maxIntron=1500 -noHead
            less -S primer.fa.psl |perl -ne '{@a=split;\$stat{\$a[9]}++ if(\$a[0]/(\$a[10]-10) > 0.95)}END{foreach \$a(keys %stat){print \$a if(\$stat{\$a}>1)}}' > multi.list
            less -S ${method}.tmp.sanger.indel.result.xls multi.list|perl -ne 'if(/Chrom/){print \$_}else{chomp;@a=split;if(scalar @a > 4){\$id=join(",",\$a[0],\$a[1],\$a[-2],\$a[-1])}else{\$id=\$_};\$stat{\$id}++;\$info{\$id}=\$_;if(\$stat{\$id} == 1){print \$info{\$id},"\\n"}}' >  ${method}.sanger.indel.result.xls
            """
    }
}

process filter_and_norm{
    publishDir "${params.outdir}/03.file_treatment/overlap_filter_norm", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    tag "${method}"
    cpus 1
    memory "8G"
    input:
        tuple val(method), path(extract_vcf_file) from extract_vcf_ch2
        file rawvcf from raw
        file ref_file from ref
    output:
        tuple val("${method}"), path("${method}.overlap_filter_norm.vcf.gz") into norm_vcf_ch
    script:
        """
        ${source_env}
        ${export_path}
        bgzip -c ${extract_vcf_file} > ${method}.split.vcf.gz
        tabix ${method}.split.vcf.gz || tabix -f -C ${method}.split.vcf.gz
        bcftools query -f "%CHROM\\t%POS\\n" ${extract_vcf_file}|perl -ne '@a=split;if(\$a[1] < 31){\$a[1] = 31;}print join("\\t",@a[0],\$a[1]-31,\$a[1]+30),"\\n";' > extract.bed
        bcftools query -T extract.bed -f "%CHROM\\t%POS\\n" ${rawvcf}| awk -v OFS="\t" '{print \$1,\$2-1,\$2}' > raw.bed
        bedtools intersect -a extract.bed -b raw.bed  -c|awk '\$4 == 1' > ${method}.overlap_filter.bed
        bcftools view -R ${method}.overlap_filter.bed ${method}.split.vcf.gz -o ${method}.overlap_filter.vcf.gz
        bcftools norm -f ${ref_file} -o ${method}.overlap_filter_norm.vcf.gz  -m - ${method}.overlap_filter.vcf.gz
        """
}

process split_vcf{
    publishDir "${params.outdir}/03.file_treatment/split_vcf", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    tag "${method}"
    input:
        tuple val(method), path(norm_vcf_file) from norm_vcf_ch
    output:
        path "*" into split_vcf_ch
    script:
        """
        ${source_env}
        gzip -d -c ${norm_vcf_file} > ${method}.norm.vcf
        java -jar ${env_path}/snpEff/SnpSift.jar split -l 1000 ${method}.norm.vcf
        rm -f ${method}.norm.vcf
        """
}

split_vcf_ch2=split_vcf_ch.flatten()


process get_flank_seq{
    publishDir "${params.outdir}/03.file_treatment/t_polymarker_input", patetrn: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    tag "${split_vcf_file}"
    input:
        each split_vcf_file from split_vcf_ch2
        file ref_file from ref
    output:
        file "${key}.flank_seq.fa"
        tuple val(method), val(key), path("${key}.for_blast.fa"), path("${key}.t_polymarker_input.csv") into t_polymarker_input_ch
    script:
        key=split_vcf_file.name.split("\\.")[0]+"_"+split_vcf_file.name.split("\\.")[2]
        method=split_vcf_file.name.split("\\.")[0]
        """
        ${source_env}
        gzip -c ${split_vcf_file} > ${key}.norm.vcf.gz
        python3 ${baseDir}/bin/get_pos.py ${key}.norm.vcf.gz ${key}.overlap_filter_norm.list
        python3 ${baseDir}/bin/get_flank_bed.py --infile ${key}.overlap_filter_norm.list --outfile ${key}.flank.bed --flank_len 50
        seqtk subseq ${ref_file} ${key}.flank.bed > ${key}.flank_seq.fa
        python3 ${baseDir}/bin/get_input_config.py --infile ${key}.flank_seq.fa --outfile ${key}.t_polymarker_input.csv --snp_list ${key}.overlap_filter_norm.list
        python3 ${baseDir}/bin/parse_polymarker_input.py ${key}.t_polymarker_input.csv ${key}.for_blast.fa
        """
}

process blast{
    publishDir "${params.outdir}/04.blast/blast_out", pattern: "*.blast_out.txt"
    publishDir "${params.outdir}/04.blast/flank_ranges", pattern: "*.flank_ranges.txt"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory "32G"
    tag "${method}_${key}"
    input:
        file ref_file from ref
        tuple val(method), val(key), path(for_blast_fa), path(t_polymarker_input) from t_polymarker_input_ch
    output:
        file "${key}.blast_out.txt"
        tuple val(method), val(key), path("${key}.final.ranges.txt") into flank_ranges_raw_ch, flank_ranges_raw_ch2
    script:
        """
        ${source_env}
        ${export_path}
        makeblastdb -in ${ref_file} -dbtype nucl -parse_seqids
        blastn -task blastn -db ${ref_file} -query ${for_blast_fa} -outfmt "6 std qseq sseq slen" -word_size 11 -num_threads 4 -out ${key}.blast_out.txt
        python3 ${baseDir}/bin/getflanking.py ${t_polymarker_input} ${key}.blast_out.txt ${key}.flank_ranges.txt
        # less -S ${key}.flank_ranges.txt | perl -ne '{chomp;@a=split;\$stat{\$a[0]}++;\$info{\$a[0]}=\$_;}END{foreach \$a(keys %stat){print \$info{\$a},"\\n" if(\$stat{\$a} == 1);}}' > ${key}.final.ranges.txt
        less -S ${key}.flank_ranges.txt > ${key}.final.ranges.txt
        """
}

flank_ranges_ch=flank_ranges_raw_ch
    .filter(v -> v[2].size()>0)

flank_ranges_ch2=flank_ranges_raw_ch2
    .filter(v -> v[2].size()>0)


process kasp{
    publishDir "${params.outdir}/05.kasp/${method}", pattern: "KASP_output/selected_KASP_primers*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory "8G"
    tag "${method}_${key}"
    input:
        file ref_file from ref
        tuple val(method), val(key), path(flank_ranges) from flank_ranges_ch
    output:
        tuple val(method), file("KASP_output/selected_KASP_primers*") into kasp_result_ch
    script:
        """
        ${source_env}
        ${export_path}
        python3 ${baseDir}/bin/get_flank_bed2.py --infile ${flank_ranges}
        gawk '{ print \$2,\$3,\$4 > "temp_marker_"\$1".txt" }' ${flank_ranges}
        for i in flank_bed/*
        do
        base_name=`basename \${i}`
        seqtk subseq ${ref_file} \${i} > flanking_temp_marker_\${base_name}.fa
        done
        python3 ${baseDir}/bin/getkasp3.py 63 25 0
        """
}


process caps{
    publishDir "${params.outdir}/06.caps/${method}", pattern: "CAPS_output/selected_CAPS_primers*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory "8G"
    tag "${method}_${key}"
    input:
        file ref_file from ref
        tuple val(method), val(key), path(flank_ranges) from flank_ranges_ch2
    output:
        tuple val(method), file("CAPS_output/selected_CAPS_primers*") into caps_result_ch
    script:
        """
        ${source_env}
        ${export_path}
        python3 ${baseDir}/bin/get_flank_bed2.py --infile ${flank_ranges}
        gawk '{ print \$2,\$3,\$4 > "temp_marker_"\$1".txt" }' ${flank_ranges}
        for i in flank_bed/*
        do
        base_name=`basename \${i}`
        seqtk subseq ${ref_file} \${i} > flanking_temp_marker_\${base_name}.fa
        done
        python3 ${baseDir}/bin/getCAPS.py 200 63 25 0
        """
}


kasp_result_ch2=kasp_result_ch
.transpose()
.groupTuple(by:0)

caps_result_ch2=caps_result_ch
.transpose()
.groupTuple(by:0)

process arrange_kasp_result{
    publishDir "${params.outdir}/05.kasp", pattern: "${method}.Potential_KASP_primers.tsv"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    tag "${method}"
    input:
        tuple val(method), file(kasp_output) from kasp_result_ch2
    output:
        tuple val(method), file("${method}.Potential_KASP_primers.tsv") into kasp_primers_ch
    script:
        """
        cat selected_*.txt > ${method}.Potential_KASP_primers.tsv
        """
}


process arrange_casp_result{
    publishDir "${params.outdir}/06.caps", pattern: "${method}.Potential_CAPS_primers.tsv"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    tag "${method}"
    input:
        tuple val(method), file(caps_output) from caps_result_ch2
    output:
        tuple val(method), file("${method}.Potential_CAPS_primers.tsv") into caps_primers_ch
    script:
        """
        cat selected_*.txt > ${method}.Potential_CAPS_primers.tsv
        """
}

kasp_merge_ch=kasp_primers_ch
.collect()
.flatten()
.collate(2)

caps_merge_ch=caps_primers_ch
.collect()
.flatten()
.collate(2)

kasp_caps_primers_ch=kasp_merge_ch
.combine(caps_merge_ch, by:0)

process result{
    publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.caps.result.xls"
    publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.dcaps.result.xls"
    publishDir "${params.outdir}/07.result/${method}", pattern: "${method}.kasp.result.xls"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "4G"
    input:
        tuple val(method), file(kasp_primers), file(caps_primers) from kasp_caps_primers_ch
    output:
        file "${method}.caps.result.xls"
        file "${method}.dcaps.result.xls"
        file "${method}.kasp.result.xls"
        """
        ${source_env}
        ${export_path}
        python3 ${baseDir}/bin/capskasp_trimlist.py ${caps_primers} ${method}.caps_primers_result.xls ${kasp_primers} ${method}.kasp.result.xls
        python3 ${baseDir}/bin/split_caps_dcaps_result.py --infile ${method}.caps_primers_result.xls --caps_result ${method}.caps.result.xls --dcaps_result ${method}.dcaps.result.xls
        """
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
