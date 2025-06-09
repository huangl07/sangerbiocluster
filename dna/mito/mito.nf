
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
            tuple val(name),read1,read2,val(type) from reads
        output:
            tuple val(name),"${name}_clean.1.fastq.gz","${name}_clean.2.fastq.gz",val(type) into clean_reads
            file("${name}.qc.stat")  into qcstat
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
            tuple val(name),read1,read2,val(type) from clean_reads
        output:
            tuple val(name),"${name}.realign.bam","${name}.realign.bam.bai" into align_bam
            file "*"
            file "${name}.metris" into metris
            file "${name}.depth" into depth
            file "${name}.chrM.coverage" into coverage
            file "${name}.flagstat" into mapstat
        script:
        if( type.grep("DC") ){
            """
            source ~/app/bioinfo/dna/new.rc
            samtools dict ${params.ref} > ref.dict
            bwa mem -M -a -t 16 -R "@RG\\tID:${name}\\tLG:${name}\\tLB:2\\tPL:illumina\\tSM:${name}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ${params.ref} ${read1} ${read2}| samtools view  -bS - > ${name}.mapping.bam
            samtools sort  -@ 16 ${name}.mapping.bam -o ${name}.sorted.bam
            samtools index ${name}.sorted.bam
            ln -s ${name}.sorted.bam ${name}.mkdup.bam
            ln -s ${name}.sorted.bam.bai ${name}.mkdup.bam.bai
            echo "0\t0\t0\t0\t0\t0\n" > ${name}.metris
            ~/app/bioinfo/dna/sentieon-genomics-202112.06/bin/sentieon driver -r ${params.ref}  -i ${name}.mkdup.bam  --algo Realigner ${name}.realign.bam
            samtools depth -b ${params.bed} ${name}.realign.bam > ${name}.depth
            samtools coverage -r chrM:0-16570 ${name}.realign.bam  > ${name}.chrM.coverage
            sed -i 's/chrM/${name}\tchrM/' ${name}.chrM.coverage
            samtools flagstats  ${name}.realign.bam  > ${name}.flagstat
            """
        }else{
            """
            source ~/app/bioinfo/dna/new.rc
            export SENTIEON_LICENSE=/mnt/lustre/users/sanger-dev/app/bioinfo/WGS/sentieon-genomics-201808/MajorBio_cluster_0.37.lic
            samtools dict ${params.ref} > ref.dict
            bwa mem -M -a -t 16 -R "@RG\\tID:${name}\\tLG:${name}\\tLB:2\\tPL:illumina\\tSM:${name}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq" ${params.ref} ${read1} ${read2}| samtools view  -bS - > ${name}.mapping.bam
            samtools sort  -@ 16 ${name}.mapping.bam -o ${name}.sorted.bam
            samtools index ${name}.sorted.bam
            java -Xmx8G -Djava.io.tmpdir=./temp/ -jar /mnt/lustre/users/sanger-dev/app/bioinfo/gene-structure/picard.jar MarkDuplicates INPUT=${name}.sorted.bam OUTPUT=${name}.mkdup.bam METRICS_FILE=${name}.metris VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 
            samtools index ${name}.mkdup.bam
            ~/app/bioinfo/dna/sentieon-genomics-202112.06/bin/sentieon driver -r ${params.ref}  -i ${name}.mkdup.bam  --algo Realigner ${name}.realign.bam
            samtools depth -b ${params.bed} ${name}.realign.bam > ${name}.depth
            samtools coverage -r chrM:0-16570 ${name}.realign.bam  > ${name}.chrM.coverage
            sed -i 's/chrM/${name}\tchrM/' ${name}.chrM.coverage
            samtools flagstats  ${name}.realign.bam  > ${name}.flagstat
            """
        }
}
process mtDNA{
    publishDir "${params.outdir}/03.mtDNA" , pattern: "*"
		tag "mtDNA"
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
        bcftools mpileup --threads 8 -d 30000 -Ou -f ${params.ref} -R ${params.bed} ${bam} -o ${name}.mpileup.bcf
        bcftools call --threads 8 -mA -Ob -o ${name}.call.bcf ${name}.mpileup.bcf
        bcftools index ${name}.call.bcf
        bcftools norm --threads 8 -f ${params.ref}  ${name}.call.bcf -Ob -o ${name}.norm.bcf
        bcftools filter --threads 8 --IndelGap 5   ${name}.norm.bcf -Ob -o ${name}.filter.bcf
        bcftools index ${name}.filter.bcf
        perl /mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/mito/bin/bcf.pl -b  ${name}.filter.bcf -o ${name}
        cat ${name}.consensus.fa ${baseDir}/bin/mtDNA.fa > align.fasta
        muscle -in align.fasta -out ${name}.align -clw
        """
}

process stat{
    publishDir "${params.outdir}/04.stat" , pattern: "*"
		tag "stat"
        cpus 16
        executor "slurm"
        queue "SANGERDEV"
        memory "10G"

        input:
            file(qc) from qcstat.collect()
            file(depth) from depth.collect()
            file(coverage) from coverage.collect()
            file(mapstat) from mapstat.collect()
            file(metris) from metris.collect()
        output:
            file "stat.xls"
        script:
        """
        perl ${baseDir}/bin/map_stat.pl -b ./ -o stat.xls
        """
}