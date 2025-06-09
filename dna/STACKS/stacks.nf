#!/usr/bin/env nextflow
params.fqlist = 'fq.list'
params.outdir = "demo"
params.help = false
params.ploid=2
params.method="GATK"
params.queue="SANGERDEV"
def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --fqlist    <file>  rawdata file list
    --outdir    <dir>   output dir
    --queue     <str>   SANGERDEV or DNA


    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

process readcsv{
    publishDir "${params.outdir}/00.prepare", pattern:"*.csv"
    executor "local"
    input:
        file csv from file(params.fqlist)
    output:
        file "newfqlist.csv" into newfqlist
    script:
    """
        less -S ${csv}|perl -ne '@a=split;\$stat{\$a[0]}=1;chomp;print \$_,"\\t",scalar keys %stat,"\\n";' > newfqlist.csv
    """
}
reads_count = newfqlist
              .splitCsv(header:false,sep:'\t')
              .groupTuple().view()

process mergefq{
    publishDir "${params.outdir}/01.mergeFQ", pattern:"*.gz"
    queue "${params.queue}"
    executor "slurm"
    input:
        tuple val(name),read1,read2,val(id) from reads_count
    output:
        tuple name,"${name}.R1.fastq.gz","${name}.R2.fastq.gz",id into reads_ch_2
        file "${name}*.fastq.gz"
    script:
        fqsize=read1.size()
        fq1=read1.join(" ")
        fq2=read2.join(" ")
    if (fqsize > 1)
        """
            cat ${fq1} > ${name}.R1.fastq.gz
            cat ${fq2} > ${name}.R2.fastq.gz
        """
    else
        """
            ln -s ${read1[0]} ${name}.R1.fastq.gz
            ln -s ${read2[0]} ${name}.R2.fastq.gz
        """
}

process cleanfq{
    publishDir "${params.outdir}/02.clean", pattern:"*.gz"
    publishDir "${params.outdir}/09.result/fig/", pattern:"*.png"
    publishDir "${params.outdir}/09.result/fig/", pattern:"*.pdf"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "200G"
    input:
        tuple val(name),read1,read2,val(id) from reads_ch_2
    output:
        tuple name,"${name}.1.fq.gz","${name}.2.fq.gz",id into reads,reads_ch_3
        file "${name}.stat" into jsons
        file "*"
    script:
    """
        #source ~/app/bioinfo/dna/new.rc
        fastp -i ${read1} -o ${name}.1.fq.gz -I ${read2} -O ${name}.2.fq.gz -w 8 -l 140 --max_len1 140 --max_len2 140 -j ${name}.json
        perl ${baseDir}/bin/fastp.pl -i ${name}.json -o ${name}
        Rscript ${baseDir}/bin/ngsqc.r --base ${name}.raw.atgcn --qual ${name}.raw.qual --key ${name}.raw --od ./
        Rscript ${baseDir}/bin/ngsqc.r --base ${name}.clean.atgcn --qual ${name}.clean.qual --key ${name}.clean --od ./
    """
}

process ustacks{
    publishDir "${params.outdir}/03.ustacks", pattern:"*.gz"
    queue "${params.queue}"
    executor "slurm"
    cpus 8
    memory "200G"
    input:
        tuple val(name),read1,read2,val(id) from reads
    output:
        tuple val(name),"${name}.tags.tsv.gz","${name}.snps.tsv.gz","${name}.alleles.tsv.gz",read1,read2 into ustacks,ustacks1,ustacks2
        file "${name}.tags.tsv.gz" into tags_ch,tags_ch1,tags_ch2
        file "${name}.snps.tsv.gz" into snps_ch,snps_ch1,snps_ch2
        file "${name}.alleles.tsv.gz" into alleles_ch,alleles_ch1,alleles_ch2
        file "*"
    script:
    """
        #source ~/app/bioinfo/dna/new.rc
        ustacks -t gzfastq -f ${read1} -o ./ --name ${name} -i ${id[0]} -p 8
    """
}

process cstacks{
    publishDir "${params.outdir}/04.cstacks",pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "1000G"
    input:
        file tags from tags_ch.collect()
        file snps from snps_ch.collect()
        file alleles from alleles_ch.collect()
    output:
        file "catalog.alleles.tsv.gz" into cstacks_alleles,cstacks_alleles1,cstacks_alleles2
        file "catalog.snps.tsv.gz" into cstacks_snps,cstacks_snps1,cstacks_snps2
        file "catalog.tags.tsv.gz" into cstacks_tags,cstacks_tags1,cstacks_tags2
        file "*"
    script:
         def stacks = tags.collect{ "$it" }.join(' ')
    """
        #source ~/app/bioinfo/dna/new.rc
        mkdir cstacks
        ls *.tags.tsv.gz|perl -ne 'chomp;\$a=\$_;\$a=~s/\\.tags\\.tsv\\.gz//g;\$n++;if(\$n == 1){print "cstacks  -o ./ -s \$a -n 4 -p 16\\n"}else{print "cstacks  -o ./ --catalog ./catalog -s \$a -n 4 -p 16 && touch \$n\\.log\\n"}' >cstacks.sh
        sh cstacks.sh
    """
}

process sstacks{
    publishDir "${params.outdir}/05.sstacks",pattern:"*.gz"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "300G"
    input:
        file calleles from cstacks_alleles
        file csnps from cstacks_snps
        file ctags from cstacks_tags
        tuple val(name),tags,snps,alles,read1,read2 from ustacks1
    output:
        tuple val(name),tags,snps,alles,"${name}.matches.tsv.gz",read1,read2 into sstacks,sstacks_ch1
        file "${name}.matches.tsv.gz" into sstacks_ch
        file "*"
    script:
    """
        ln -s ${tags}
        ln -s ${snps}
        ln -s ${alles}
        ln -s ${read1}
        ln -s ${read2}
        sstacks -c ./ -s ${name} -o ./ -p 8
    """
}
process tsv2bam{
    publishDir "${params.outdir}/06.tsv2bam",pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "200G"
    input:
        tuple val(name),tags,snps,alles,matches,read1,read2 from sstacks
        file calleles from cstacks_alleles1
        file csnps from cstacks_snps1
        file ctags from cstacks_tags1
    output:
        tuple name,"${name}.matches.bam","${name}.depth" into tsv2bam
        file "${name}.matches.bam" into bams
        file "${name}.depth" into depth
        file "*"
    script:
    """
        ln -s ${tags}
        ln -s ${snps}
        ln -s ${alles}
        ln -s ${matches}
        ln -s ${read1}
        ln -s ${read2}
        tsv2bam -P ./ -s ${name} -t 8
        samtools coverage ${name}.matches.bam > ${name}.depth
    """
}
process genotypes{
    publishDir "${params.outdir}/07.populations",pattern:"*"
    publishDir "${params.outdir}/09.results",pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    cpus 16
    memory "200G"
    input:
        tuple val(name),bams,depths from tsv2bam.collect()
        file tags from tags_ch1.collect()
        file snps from snps_ch1.collect()
        file alleles from alleles_ch1.collect()
        file calleles from cstacks_alleles2
        file csnps from cstacks_snps2
        file ctags from cstacks_tags2
        file sstacks from sstacks_ch.collect()
        file bams from bams.collect()
        file jsons from jsons.collect()
        file depth from depth.collect()
    output:
        file "populations.final.vcf.gz" into vcf
        file "populations.final.vcf.gz"
        file "populations.final.vcf.gz.tbi"
        file "sample.stat"
        file "snp.stat.xls"
        file "populations.loci.fa"
        file "pop.qc"
        file "snp.stat.all.xls"
        file "populations.loci.fa.fai"
    script:
    """
        ls *.matches.bam|perl -ne 'chomp;s/\\.matches\\.bam//g;print \$_,"\\t","pop1\\n";' > uniq.txt
        cat *.stat|sort|uniq > pop.qc
        gstacks -P ./ -O ./ -t 8 -M uniq.txt
        populations -P ./ -M uniq.txt -t 16 -O ./ --vcf --fasta-loci
        less -S populations.snps.vcf|perl -ne 'if(/^#/){print \$_}else{\$n++;@a=split;print "sca1\\t\$n\\t",join("\\t",@a[2..\$#a]),"\\n"}' > populations.final.vcf
        snpnum=`grep -v "#" populations.final.vcf|wc -l`
        sed -i '14 i ##contig=<ID=sca1,length='\${snpnum}',assembly=unknown>' populations.final.vcf
        bgzip -@ 8 populations.final.vcf
        tabix populations.final.vcf.gz
        perl ${baseDir}/bin/stacks.pl -d ./ -o ./
        bcftools +smpl-stats populations.final.vcf.gz -o pop.stat
        echo -e 'sample\\tnum_gt\\tnum_nonref_gt\\tnum_homo_ref\\tnum_homo_alt\\tnum_het\\tnum_snv\\tnum_indel\\tnum_singleton\\tnum_miss\\tnum_ts\\tnum_tv\\tts/tv' > snp.stat.xls
        grep 'FLT0' pop.stat | sed '1d' |cut -f 2-7,9- >> snp.stat.xls
        echo -e 'total\\tsnvs\\tindels\\tsingletons\\ttransitions\\ttransversions\\tts/tv' > snp.stat.all.xls
        grep 'SITE0' pop.stat | cut -f 2- >> snp.stat.all.xls
        sed -i "1d" populations.loci.fa
        samtools faidx populations.loci.fa
    """
}
