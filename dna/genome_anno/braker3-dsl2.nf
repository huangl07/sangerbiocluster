#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

params.ref="ref.fa"
params.db="/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Arthropoda.fa"
params.fqlist="None"
params.outdir="result"
params.help = false

def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run pasa.nf --ref genome --outdir ./result

    --outdir    <dir>   output dir
    --ref   <file>  reference genome fasta file
    --fqlist   <file>  sample\\tfq1\\tfq2，在提供fq的情况下，运行Braker3_bam
    --db    <file>  database, 目前提供5个
    1、节肢动物门：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Arthropoda.fa
    2、脊椎动物门：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Vertebrata.fa
    3、植物界：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Viridiplantae.fa
    4、真菌：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Fungi.fa
    5、真核生物：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Eukaryota.fa
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


process RepeatModeler {
	    publishDir "${params.outdir}/01.RepeatAnno", pattern: "*"
		tag "RepeatModeler"
        queue "SANGERDEV"
        cpus 32
        memory "400G"
        executor "slurm"
        input:
            path ref
        output:
            path "genome-families.fa" , emit: ref_mode

        script:
        """
            source ~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
            conda activate ~/app/bioinfo/dna/miniconda3/envs/repeat
            BuildDatabase -name genome ${ref}
            RepeatModeler -database genome  -threads 32 -LTRStruct
        """
}

process RepeatMasker {
	    publishDir "${params.outdir}/01.RepeatAnno", pattern: "*"
		tag "RepeatMasker"
        queue "SANGERDEV"
        cpus 40
        memory "400G"
        executor "slurm"
        input:
            path ref
            path mode
        output:
            path "*.masked" , emit: ref_masked
            path "*.out" , emit: ref_masked_table
            path "*.tbl" , emit: ref_masked_stat

        script:
        """
            source ~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
            conda activate ~/app/bioinfo/dna/miniconda3/envs/repeat
            RepeatMasker -pa 40 -lib ${mode} ${ref} -small -xsmall -gff -html -poly
        """
}


process Braker3_only {
        publishDir "${params.outdir}/02.braker3", pattern: "*"
		tag "braker3_nobam"
        queue "SANGERDEV"
        cpus 40
        memory "400G"
        executor "slurm"
        input:
            path ref
            path db
        output:
            path "braker.*"

        script:
        """
            source ~/app/bioinfo/dna/new.rc
            set -h
            source ~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
            conda activate ~/app/bioinfo/dna/miniconda3/envs/braker3-2023-12-11/
            export GENEMARK_PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/GeneMark-ETP/bin
            export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/GeneMark-ETP/tools:\$PATH
            export PROTHINT_PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/ProtHint/bin
            less -S ${ref}|perl -ne 'chomp;@a=split;print \$a[0],"\\n";' > ref.nospace.fa
            braker.pl --genome=ref.nospace.fa --prot_seq=${db} --workingdir=result --threads 40
            ln result/braker.aa ./braker.aa
            ln result/braker.codingseq ./braker.codingseq
            ln result/braker.gtf ./braker.gtf
        """
}

process Hisat2_index {
        publishDir "${params.outdir}/02.braker3", pattern: "*"
		tag "Hisat2_index"
        queue "SANGERDEV"
        cpus 16
        memory "400G"
        executor "slurm"
        input:
            path ref
        output:
            path "genome.*", emit: genome
        script:
        """
            source ~/app/bioinfo/dna/new.rc
            ~/app/bioinfo/dna/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2-build -p 16 ${ref} genome
        """
}

process Hisat2 {
        publishDir "${params.outdir}/02.braker3", pattern: "*"
		tag "Hisat2_${name}"
        queue "SANGERDEV"
        cpus 16
        memory "120G"
        executor "slurm"
        input:
            tuple val(name),path(fq1),path(fq2)
            path genome
        output:
            path "${name}.bam*", emit: bam
        script:
        """
            source ~/app/bioinfo/dna/new.rc
            ~/app/bioinfo/dna/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -x genome -1 ${fq1} -2 ${fq2} -S ${name}.sam --rna-strandness RF -p 16 --dta
            samtools sort -@ 10 -m 2G -o ${name}.bam ${name}.sam
            rm ${name}.sam
            samtools index -@ 10 ${name}.bam || samtools index -c -@ 10 ${name}.bam
        """
}

process Braker3_bam {
        publishDir "${params.outdir}/02.braker3", pattern: "*"
		tag "braker3_bam"
        queue "SANGERDEV"
        cpus 40
        memory "400G"
        executor "slurm"
        input:
            path ref
            path db
            path "*"
        output:
            path "braker.*"

        script:
        """
            source ~/app/bioinfo/dna/new.rc
            set -h
            source /mnt/lustre/users/sanger-dev/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
            conda activate ~/app/bioinfo/dna/miniconda3/envs/braker3-2023-12-11/
            export GENEMARK_PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/GeneMark-ETP/bin
            export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/GeneMark-ETP/tools:\$PATH
            export PROTHINT_PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/ProtHint/bin
            less -S ${ref}|perl -ne 'chomp;@a=split;print \$a[0],"\\n";' > ref.nospace.fa
            bamlist=\$(ls *.bam | tr "\\n" "," | sed 's/,\$//g')
            braker.pl --genome=ref.nospace.fa --prot_seq=${db}  --workingdir=result --bam=\${bamlist} --threads=40
            ln result/braker.aa ./braker.aa
            ln result/braker.codingseq ./braker.codingseq
            awk -v FS="\\t" -v OFS="\\t" '{if (\$3=="gene") {\$9="gene_id \\"\\"\$9\\"\";"} else if (\$3=="transcript") {tmp=\$9; \$9="transcript_id \\"\\"tmp\\"\"; gene_id \\"\\"gensub(/.t.+/,\"\",\"g\",tmp)\\"\";"} print \$0}' result/braker.gtf > ./braker.gtf
        """
}

workflow {
    RepeatModeler(file(params.ref))
    RepeatMasker(file(params.ref), RepeatModeler.out.ref_mode)
    if(params.fqlist == "None"){
        Braker3_only(RepeatMasker.out.ref_masked, file(params.db))
    }
    else{
        fqch = Channel.from(file(params.fqlist))
              .splitCsv(header:false,sep:'\t')
        Hisat2_index(RepeatMasker.out.ref_masked)
        Hisat2(fqch, Hisat2_index.out.genome)
        Braker3_bam(RepeatMasker.out.ref_masked, file(params.db), Hisat2.out.bam.flatten().collect())
    }
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
