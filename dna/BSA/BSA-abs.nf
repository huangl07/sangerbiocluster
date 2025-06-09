#!/usr/bin/env nextflow
params.outdir = "demo"
params.winsize="1000000"
params.stepsize="10000"
params.help = false
params.popt="F2"
params.bulksize=30
params.bootstrap=1000
params.pvalue = 0.001
params.quantile = 0.999
params.grade=false
params.mutmap=false
params.csi=false
params.nodeepbsa=false
params.noparent=false
params.abs=false
params.env="~/app/bioinfo/dna/new.rc_1"
params.queue="SANGERDEV"
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
group=file(params.group)
chr=file(params.chr)
gff_file=file(params.gff)
anno_summary_file=file(params.anno_summary)

if(!params.grade){
process vcf2table{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory '32GB'
    input:
        file vcf from vcf
        file chr from chr
        file group from group
        file anno from anno_summary_file
    output:
        file "pop.index" into index1,index2,index3,index4,index5
        file "pop.table" into table_file, pop_table_ch2
        file "pop.final.anno" into anno_table
        file "*"
    script:
        """
        source ${params.env}
        perl ${baseDir}/bin/vcf2table.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
        perl ${baseDir}/bin/bsa_calc.pl --table pop.table --group ${group} --popt ${params.popt} -out pop.index
        perl ${baseDir}/bin/anno.pl --table pop.anno --anno ${anno} -out pop.final.anno
        python3 ${baseDir}/bin/bsa_stat_filtered_pop_table.py --infile pop.table --outfile snp_indel_gene.stat --chrfile ${chr}
        """
}

if(!params.csi){
process abstract_filter{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
        file pop_table from pop_table_ch2
    output:
        file "filter.vcf.gz" into pop_filtered_vcf_ch3
        file "filter.vcf" into pop_filtered_vcf_ch, pop_filtered_vcf_ch2
    script:
    """
        source ${params.env}
        cut -f 1,2 ${pop_table} | sed '1d' > pos_list.txt
        tabix ${vcf}
        python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list
        bcftools view -S mix.list -R pos_list.txt ${vcf} -o filter.vcf -Ov --threads 8
        bgzip -k filter.vcf
    """
}
}else{
process abstract_filter_csi{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
        file pop_table from pop_table_ch2
    output:
        file "filter.vcf.gz" into pop_filtered_vcf_ch3
        file "filter.vcf" into pop_filtered_vcf_ch, pop_filtered_vcf_ch2
    script:
    """
        source ${params.env}
        cut -f 1,2 ${pop_table} | sed '1d' > pos_list.txt
        tabix -C ${vcf}
        python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list
        bcftools view -S mix.list -R pos_list.txt ${vcf} -o filter.vcf -Ov --threads 8
        bgzip -k filter.vcf
    """
}
}

process indexslid{
    publishDir "${params.outdir}/02.index-slid", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file index_file from index1
    output:
        tuple val("index"),file("pop.index.region") into index_region
        file "*"
    script:
    if(params.mutmap){
      """
        source ${params.env}
        Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue} --mutmap
        Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10
        echo index_winsize:${params.winsize} > index.params.log
        echo index_stepsize:${params.stepsize} >> index.params.log
        echo index_bootstrap:${params.bootstrap} >> index.params.log
        echo index_p:${params.pvalue} >> index.params.log
      """
    }else if(params.noparent){
        """
        source ${params.env}
        Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp --abs yes
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10
        echo index_winsize:${params.winsize} > index.params.log
        echo index_stepsize:${params.stepsize} >> index.params.log
        echo index_bootstrap:${params.bootstrap} >> index.params.log
        echo index_p:${params.pvalue} >> index.params.log
        """
    }else if(params.abs){
      """
        source ${params.env}
        Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp --abs yes
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10
        echo index_winsize:${params.winsize} > index.params.log
        echo index_stepsize:${params.stepsize} >> index.params.log
        echo index_bootstrap:${params.bootstrap} >> index.params.log
        echo index_p:${params.pvalue} >> index.params.log
      """
    }else {
        """
        source ${params.env}
        Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index"
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number 10
        echo index_winsize:${params.winsize} > index.params.log
        echo index_stepsize:${params.stepsize} >> index.params.log
        echo index_bootstrap:${params.bootstrap} >> index.params.log
        echo index_p:${params.pvalue} >> index.params.log
        """
    }
}
process loesscalc{
    publishDir "${params.outdir}/05.loess", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file table_file from index2
    output:
        file "*"
        tuple val("loess"),file("pop.loess.region") into loess_region
    script:
    """
        source ${params.env}
        Rscript ${baseDir}/bin/loess-index-snow.R --infile ${table_file} --out pop.index.loess
        Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.index.loess --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue}
        Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result pop.bootstrap.result --pcol delta --lcol slidingD --xlab chromosome --ylab "loess" --threshold ${params.quantile}  --output loess --chr ${chr} --mutmap FALSE
        Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --quantile ${params.quantile} --outfile pop.loess.region --number 10
        echo loess_q:${params.quantile} > loess.params.log
    """
}
    if(!params.mutmap){
        process Gcalc{
            publishDir "${params.outdir}/03.Gprime", pattern:"*"
            queue "${params.queue}"
            executor "slurm"
            input:
                file index_file from index3
            output:
                file "*"
                tuple val("Gprime"),file("pop.Gprime.region") into Gprime_region
            script:
            """
                source ${params.env}
                Rscript ${baseDir}/bin/G-calc.R --infile ${index_file} --outfile pop.Gprime --winsize ${params.winsize} --method bp --threshold ${params.pvalue}
                Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result pop.Gprime --pcol G --lcol Gprime --threshold ${params.quantile} --xlab chromosome --ylab "Gprime" --output Gprime  --chr ${chr}
                Rscript ${baseDir}/bin/region.R --infile pop.Gprime --ccol X.chr --pos pos --loess Gprime --quantile ${params.quantile} --outfile pop.Gprime.region --number 10
                echo gprime_winsize:${params.winsize} > gprime.params.log
                echo gprime_q:${params.quantile} >> gprime.params.log
            """
        }
        process EDslid{
        publishDir "${params.outdir}/04.ED-slid", pattern:"*"
        queue "${params.queue}"
        executor "slurm"
        input:
            file index_file from index4
        output:
            file "*"
            tuple val("ED"),file("pop.ED.region") into ED_region
        script:
        """
            source ${params.env}
            Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method bp
            Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.sliding.result --chr ${chr} --output pop.ED  --pcol slidingED --threshold ${params.quantile} --xlab chromosome --ylab "ED"
            Rscript ${baseDir}/bin/region.R --infile pop.sliding.detail --ccol X.chr --pos pos --loess slidingED --quantile ${params.quantile} --outfile pop.ED.region --number 10
            echo ed_winsize:${params.winsize} > ed.params.log
            echo ed_stepsize:${params.stepsize} >> ed.params.log
            echo ed_q:${params.quantile} >> ed.params.log
        """
        }
        process QTGseq{
        publishDir "${params.outdir}/07.QTG-seq", pattern:"*"
        queue "${params.queue}"
        executor "slurm"
        input:
            file index_file from index5
        output:
            file "*"
            tuple val("QTG"),file("pop.QTG.region") into QTG_region
        script:
        """
        source ${params.env}
        Rscript ${baseDir}/bin/smoothLOD.R  --infile ${index_file} --outdir ./ --winsize 100000
        Rscript ${baseDir}/bin/region.R --infile significant.result.txt --ccol X.chr --pos pos --loess SlidingD --CI CI --outfile pop.QTG.region --number 1
        mv significant.result.txt pop.QTG.result
        """
        }
        if(!params.nodeepbsa){
            process deepBSA_DL{
            publishDir "${params.outdir}/08.DeepBSA", pattern:"*"
            queue "AVX512"
            executor "slurm"
            input:
                file deepbsa_vcf from pop_filtered_vcf_ch
            output:
                file "*"
                tuple val("DeepBSA_DL"),file("pop.deepBSA_DL.region") into deepBSA_DL_region
            script:
            """
            source ${params.env}
            ln -s ~/app/bioinfo/dna/DeepBSA_linux_v1.4/bin/Models ./
            export R_HOME=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/lib64/R/
            export LD_LIBRARY_PATH=\$R_HOME/lib:\$LD_LIBRARY_PATH
            python3 ~/app/bioinfo/dna/DeepBSA_linux_v1.4/bin/main.py --i ${deepbsa_vcf} --m DL
            python3 ${baseDir}/bin/process_deepBSA.py --infile ./Results/filter/DL_values.txt --outfile DL_values.txt
            Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result DL_values.txt --chr ${chr} --output pop.deepBSA_DL --pcol value --lcol loess --threshold ${params.quantile} --xlab chromosome --ylab "deepBSA_DL"
            Rscript ${baseDir}/bin/region.R --infile DL_values.txt --ccol X.chr --pos pos --loess loess --quantile ${params.quantile} --outfile pop.deepBSA_DL.region --number 10
            """
            }
            process deepBSA_K{
            publishDir "${params.outdir}/08.DeepBSA", pattern:"*"
            queue "AVX512"
            executor "slurm"
            input:
                file deepbsa_vcf from pop_filtered_vcf_ch2
            output:
                file "*"
                tuple val("DeepBSA_K"),file("pop.deepBSA_K.region") into deepBSA_K_region
            script:
            """
            source ${params.env}
            ln -s ~/app/bioinfo/dna/DeepBSA_linux_v1.4/bin/Models ./
            export R_HOME=/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/lib64/R/
            export LD_LIBRARY_PATH=\$R_HOME/lib:\$LD_LIBRARY_PATH
            python3 ~/app/bioinfo/dna/DeepBSA_linux_v1.4/bin/main.py --i ${deepbsa_vcf} --m K --p 0
            python3 ${baseDir}/bin/process_deepBSA.py --infile ./Results/filter/K_values.txt --outfile K_values.txt
            Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result K_values.txt --chr ${chr} --output pop.deepBSA_K --pcol value --lcol loess --threshold ${params.quantile} --xlab chromosome --ylab "deepBSA_K"
            Rscript ${baseDir}/bin/region.R --infile K_values.txt --ccol X.chr --pos pos --loess loess --quantile ${params.quantile} --outfile pop.deepBSA_K.region --number 10
            """
            }
        }
    }


    Channel.from(index_region).mix(loess_region).set{regions}
    if(!params.mutmap){
        if(!params.nodeepbsa){
            Channel.from(regions).mix(ED_region,Gprime_region,loess_region,QTG_region,deepBSA_DL_region,deepBSA_K_region).set{regions}
        }else{
            Channel.from(regions).mix(ED_region,Gprime_region,loess_region,QTG_region).set{regions}
        }
    }
    Channel.from(index_region).mix(loess_region).set{regions2}
    if(!params.mutmap){
        if(!params.nodeepbsa){
            Channel.from(regions2).mix(ED_region,Gprime_region,loess_region,QTG_region,deepBSA_DL_region,deepBSA_K_region).set{regions2}
        }else{
            Channel.from(regions2).mix(ED_region,Gprime_region,loess_region,QTG_region).set{regions2}
        }
    }



}else{
process params.grade{
    publishDir "${params.outdir}/02.ridit", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
        file anno from anno_summary_file
    output:
        file "*"
        tuple val("ridit"),"pop.ridit.region" into regions
        file "pop.table" into table_file, pop_table_ch2
        file "pop.final.anno" into anno_table
    script:
    """
        source ${params.env}
        perl ${baseDir}/bin/vcf2table_multi.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
        perl ${baseDir}/bin/anno.pl --table pop.anno --anno ${anno} -out pop.final.anno
        python3 ${baseDir}/bin/bsa_stat_filtered_pop_table.py --infile pop.table --outfile snp_indel_gene.stat
        Rscript ${baseDir}/bin/ridit-noloess.R --index pop.table --out pop --pvalue ${params.pvalue} --group ${params.group} 
        Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.denoise.result --pcol logP --qcol CI --xlab chromosome --ylab "ridit" --output loess --chr ${chr}
        Rscript ${baseDir}/bin/region-ridit.R --infile pop.denoise.result --ccol X.chr --pos pos --loess logP --CI CI --outfile pop.ridit.region --number 10 --winsize ${params.winsize}
    """
}

if(!params.csi){
process abstract_filter_grade{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
        file pop_table from pop_table_ch2
    output:
        file "DeepBSA.vcf" into pop_filtered_vcf_ch, pop_filtered_vcf_ch2
    script:
    """
        source ${params.env}
        cut -f 1,2 ${pop_table} | sed '1d' > pos_list.txt
        tabix ${vcf}
        python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list
        bcftools view -S mix.list -R pos_list.txt ${vcf} -o deepBSA.vcf -Ov --threads 8
    """
}
}else{
process abstract_filter_grade_csi{
    publishDir "${params.outdir}/01.vcf2table", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        file vcf from vcf
        file group from group
        file pop_table from pop_table_ch2
    output:
        file "DeepBSA.vcf" into pop_filtered_vcf_ch, pop_filtered_vcf_ch2
    script:
    """
        source ${params.env}
        cut -f 1,2 ${pop_table} | sed '1d' > pos_list.txt
        tabix -C ${vcf}
        python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list
        bcftools view -S mix.list -R pos_list.txt ${vcf} -o deepBSA.vcf -Ov --threads 8
    """
}
}
}

process Orgdb {
    publishDir "${params.outdir}/00.orgdb", pattern:"*"
    executor 'slurm'
    queue "${params.queue}"
    cpus 2
    cache 'lenient'
    input:
        file anno_summary from anno_summary_file
    output:
        tuple path("TERM2GENE.txt"),path("TERM2NAME.txt") into term_ch
        file "genes_all.list" into all_gene_ch
        file "KEGG_rename.csv" into kegg_anno_ch
        file "GO_anno.csv" into go_anno_ch
        path(LIB) into EnrichDb
    script:
        """
        source ${params.env}
        mkdir LIB
        python3 ${baseDir}/bin/bsa_abstract_anno.py --infile ${anno_summary} --outfile go_kegg.list --go_anno GO_anno.csv --kegg_anno KEGG_rename.csv
        python3 ${baseDir}/bin/bsa_abstract_all_transcripts_from_anno_summary.py --infile ${anno_summary} --outfile genes_all.list
        Rscript ${baseDir}/bin/make_new_orgdb.R --annofile go_kegg.list --term2gene TERM2GENE.txt --term2name TERM2NAME.txt --outpath .
        R CMD INSTALL org.Ddemo.eg.db --library=./LIB/
        """
}

process enrich{
    publishDir "${params.outdir}/06.enrich", pattern:"*"
    queue "${params.queue}"
    executor "slurm"
    input:
        tuple val(method),region from regions
        tuple path(term2gene),path(term2name) from term_ch
        file all_genes_list from all_gene_ch
        path LIB from EnrichDb
        file KEGG_anno from kegg_anno_ch
        file GO_anno from go_anno_ch
        file anno from anno_table
        file anno_summary from anno_summary_file
		file vcf_file from vcf
    output:
        file "*"
    script:
    if(region.countLines() != 0){
    """
    #!/usr/bin/sh
        source ${params.env}
        tabix ${vcf_file}
        mkdir ${method}
        Rscript ${baseDir}/bin/abstract_pop_genes.R --regionfile ${region} --genefile ${all_genes_list} --outfile ${method}.genes_abstract.list
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$5}' > ${method}.degfile
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$6}' > ${method}.transcript.degfile
        cut -f5,6 ${method}.genes_abstract.list  > geneID_transcriptID
        csvtk cut -f GeneName,GeneID -d $'\t' -T pop.final.anno >geneNAME_geneID
        csvtk merge -L --na -- -Tt -f "gene_id;GeneID" geneID_transcriptID geneNAME_geneID >geneID_transcriptID_geneNAME 
        Rscript ${baseDir}/bin/enrich.R --degfile ${method}.transcript.degfile --term2genefile ${term2gene} --term2namefile ${term2name} --outname ${method} --db ${LIB} --outdir ${method}/ --title ${method} --genename geneID_transcriptID_geneNAME
        perl ${baseDir}/bin/anno_region.pl --table ${anno} --region ${region} --out ${method}/${method}.all.table
        python3 ${baseDir}/bin/bsa_stat_region.py --infile ${method}/${method}.all.table --abstract_gene ${method}.genes_abstract.list --anno_summary ${anno_summary} --outfile ${method}/${method}.region_stat.xls
		bcftools view --threads 4 -R ${region} ${vcf_file} -o ${method}.extract.raw.vcf.gz
    """
    }else{
    """
    echo "no-region,please check!" > ${method}.error.log
    """
    }
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

'''
process extract_region_from_vcf{
    publishDir "${params.outdir}/07.extract_vcf", pattern: "extract.vcf.gz"
    queue "${params.queue}"
    executor "slurm"
    cpus 4
    memory "32G"
    input:
        tuple val(method),region from regions2
        file vcf_file from pop_filtered_vcf_ch3
    output:
        tuple val(method), file("${method}.extract.vcf") into extract_vcf_ch, extract_vcf_ch2
        tuple val(method), file("${method}.extract.vcf.gz") into extract_gz_vcf_ch
    script:
        """
        source ${params.env}
        tabix -f ${vcf_file}
        bcftools view --threads 4 -R ${region} ${vcf_file} -o ${method}.extract.vcf.gz
        gunzip -k ${method}.extract.vcf.gz
        """
}

process primer3{
    publishDir "${params.outdir}/08.primer3", pattern: "*"
    queue "${params.queue}"
    executor "slurm"
    cpus 1
    memory "8G"
    input:
        tuple val(method), path(extract_vcf_file) from extract_vcf_ch
        file ref_file from ref
    output:
        path "${method}"    
    script:
        """
        source ${params.env}
        mkdir ${method}
        export PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/bin:\$PATH
        export LD_LIBRARY_PATH=/mnt/lustre/users/sanger-dev/app/gcc/5.1.0/lib64:\$LD_LIBRARY_PATH
        workdir=\$PWD
        perl ${baseDir}/bin/primer_design/1.p3in.pl -d ${extract_vcf_file} -r 200-300 -T1 57.0 -T2 63.0 -p 3 -ref ${ref_file} -o ./${method}
        cd /mnt/lustre/users/sanger-dev/app/bioinfo/WGS/primer3/src
        ./primer3_core --output \$workdir/${method}/variation.p3out \$workdir/${method}/variation.p3in
        cd \$workdir
        perl ${baseDir}/bin/primer_design/3.merge_p3out.v2.pl \$workdir/${method}/variation.p3out 3
        """
}
'''