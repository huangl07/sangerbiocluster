#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

params.ref="ref.fa"
params.transcript="Trinity.fasta"
params.outdir="result"
params.help = false

def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow-1 run -bg pasa.nf --ref genome.fa --transcript trinity.fa --outdir ./result

    --outdir    <dir>   output dir
    --ref   <file>  reference genome fasta file
    --transcript    <file>  illumina trinity fasta

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}




process Pasa {
	    publishDir "${params.outdir}", pattern: "*", mode: 'copy'
		tag "Pasa"
        queue "SANGERDEV"
        cpus 32
        memory "32G"
        executor "slurm"
        
        input:
            path ref 
            path cds 
        output:
            path "result.*"

        script:
        """
            mkdir -p temp
            mkdir -p data

            cp  ${ref} data/genome_sample.fasta
            cp  ${cds} data/all_transcripts.fasta
            cp  /mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/alignAssembly.config data/alignAssembly.config

            docker load -i /mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/pasa.image
            workpath=`pwd`
            docker run -i --rm --privileged -v \${workpath}/temp:/tmp -v \${workpath}/data:/data 87d49c3ed527 bash -c 'cd /data && /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl  -c alignAssembly.config -C -R --ALIGNER gmap -g genome_sample.fasta -t all_transcripts.fasta --CPU 32'


            cp data/sample_mydb_pasa.assemblies.fasta result.fa
            cp data/sample_mydb_pasa.pasa_assemblies.gtf result.gtf
            cp data/sample_mydb_pasa.pasa_assemblies.gff3 result.gff3
        """
}

workflow {
    Pasa(params.ref, params.transcript)
}
