
params.ref="/mnt/lustre/users/sanger-dev/ysource/files/mtDNA/files/GRCh38/Homo_sapiens_assembly38.fasta"
params.bed="/mnt/lustre/users/sanger-dev/ysource/scripts/mtDNA/mtDNA/mtDNA.bed"
params.outdir="demo"
reads = Channel.from(file(params.fqlist))
              .splitCsv(header:false,sep:'\t')
              .groupTuple()
process fastps{
    publishDir "${params.outdir}/01.cleanFQ" , pattern: "*"
		tag "fastp"
        cpus 16
        memory "10G"
        executor "slurm"
        queue "SANGERDEV"
        input:
            tuple val(name), read1,read2 from reads
        output:
            tuple val(name),"${name}_clean.1.fastq.gz","${name}_clean.2.fastq.gz" into clean_reads
            file("${name}.stat")  into qcstat
            file "${name}.json"
            file "${name}.html"
        script:
        fq1=read1.join(" ")
        fq2=read2.join(" ")
        """
            source ~/app/bioinfo/dna/new.rc
            fastp -i ${fq1} -o ${name}_clean.1.fastq.gz -I ${fq2} -O ${name}_clean.2.fastq.gz -w 8 -h ${name}.html -j ${name}.json  --cut_front --cut_tail --average_qual 35 --cut_mean_quality 30 --low_complexity_filter --qualified_quality_phred 35
            perl ${baseDir}/bin/fastp.pl -i ${name}.json -o ${name}.qc
        """
}
process mapping{
    publishDir "${params.outdir}/02.mapping" , pattern: "*"
		tag "mapping"
        cpus 16
        executor "slurm"
        queue "SANGERDEV"
        memory "20G"
        input:
            tuple val(name),read1,read2 from clean_reads
        output:
            tuple val(name),"${name}.sorted.bam","${name}.sorted.bam.bai" into align_bam
            file "*"
            file "${name}.depth" into depth
            file "${name}.chrM.coverage" into coverage
        script:
        """
        source ~/app/bioinfo/dna/new.rc
        samtools dict ${params.ref} > ref.dict
        bwa mem -M -a -t 16 -R "@RG\\tID:${name}\\tLG:${name}\\tLB:2\\tPL:illumina\\tSM:${name}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ${params.ref} ${read1} ${read2}| samtools view  -bS - > ${name}.mapping.bam
        samtools sort  -@ 16 ${name}.mapping.bam -o ${name}.sorted.bam
        samtools index ${name}.sorted.bam
        java -Xmx8G -Djava.io.tmpdir=./temp/ -jar /mnt/lustre/users/sanger-dev/app/bioinfo/gene-structure/picard.jar MarkDuplicates INPUT=${name}.sorted.bam OUTPUT=${name}.mkdup.bam METRICS_FILE=${name}.metris VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 
        samtools depth -b ${params.bed} ${name}.mkdup.bam > ${name}.depth
        samtools coverage -r chrM:1-16570 ZB_7.mkdup.bam  > ${name}.chrM.coverage
        """
}
process mtDNA{
    publishDir "${params.outdir}/03.mtDNA" , pattern: "*"
		tag "mapping"
        cpus 16
        executor "slurm"
        queue "SANGERDEV"
        memory "10G"

        input:
            tuple val(name),bam,bai from align_bam
        output:
            file "*"
        script:
        """
        source ~/app/bioinfo/dna/new.rc
        bcftools mpileup --threads 8 -d 30000 -Ou -f ${params.ref} -R ${params.bed} ${bam} -o ${name}.mpileup.bcf
        bcftools call --threads 8 -mA -Ob -o ${name}.call.bcf ${name}.mpileup.bcf
        bcftools index ${name}.call.bcf
        bcftools norm --threads 8 -f ${params.ref}  ${name}.call.bcf -Ob -o ${name}.norm.bcf
        bcftools filter --threads 8 --IndelGap 5   ${name}.norm.bcf -Ob -o ${name}.filter.vcf.gz
        tabix ${name}.filter.vcf.gz
        cat ${params.ref}|bcftools consensus ${name}.filter.vcf.gz > ${name}.consensus.fa
        bcftools call --threads 8 -mv -Ob -o ${name}.variant.bcf ${name}.mpileup.bcf
        bcftools norm --threads 8 -f ${params.ref}  ${name}.variant.bcf -Ob -o ${name}.vnorm.bcf
        bcftools filter --threads 8 --IndelGap 5   ${name}.vnorm.bcf -Ob -o ${name}.vfilter.bcf
        bcftools query -Hf "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%TGT\\t]\\t%DP\\n" ${name}.vfilter.bcf >${name}.variant.table
        seqtk subseq ${name}.consensus.fa ${params.bed} > chrm.fa
        """
}

process stat{
    publishDir "${params.outdir}/04.stat" , pattern: "*"
		tag "mapping"
        cpus 16
        executor "slurm"
        queue "SANGERDEV"
        memory "10G"

        input:
            file(qc) from qcstat.collect()
            file(depth) from depth.collect()
            file(coverage) from coverage.collect()
        output:
            file "*"
        script:
        """
        cat *.qc.stat > qc.stat
        cat *.coverage > 
        """
}