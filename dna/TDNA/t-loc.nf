#!/usr/bin/env nextflow

/*
 * Enable DSL 1 syntax
 */
nextflow.enable.dsl = 1

/*
 * Define the default parameters
 */

params.ref="ref.fa"
params.TDNA=""
params.fq1=""
params.fq2=""
params.outdir="result"
params.help = false

def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run -bg tloc.nf --ref genome --TDNA --fq1 --fq2 --outdir ./result

    --outdir    <dir>   output dir
    --ref   <file>  reference genome fasta file

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

ref_file = file(params.ref)
TDNA = file(paramas.TDNA)
fq1 = file(params.fq1)
fq2 = file(params.fq2)

process index {
	    publishDir "${params.outdir}/01.index", pattern: "*", mode: 'link'
		tag "index"
        queue "SANGERDEV"
        cpus 8
        memory "32G"
        executor "slurm"
        input:
            path ref from ref_file
        output:
            path "ref.fa*" 
            path "ref.fa*" into index
        script:
        """
            source ~/app/bioinfo/dna/new.rc
            ln -s ${ref} ref.fa
            ~/app/bioinfo/dna/env/bin/bwa-mem2 index re.fa
        """
}

process TLOC {
	    publishDir "${params.outdir}/02.t-loc", pattern: "*", mode: 'link'
		tag "TLOC"
        queue "SANGERDEV"
        cpus 8
        memory "30G"
        executor "slurm"
        input:
            path index from index
            path TDNA from TDNA
            path fq1 from fq1
            path fq2 from fq2
        output:
            path "Result/*"
            
        script:
        """
            source ~/app/bioinfo/dna/new.rc
            python2.7 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/demo/output/test_py/T-LOC-mem/T-LOC2.py  --fastq ${fq1},${fq2} --genome_Bwa ${index}  --TDNA ${TDNA} --output ./
        """
}
