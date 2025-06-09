#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.env = "~/app/bioinfo/dna/new.rc"
params.queue = "SANGERDEV"
def helpMessage() {
   
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --mark    <file>  input mark file
    --map     <file>   output dir
    --trt     <file>  input group file    
    --outdir  <dir>   output dir
    --popt  <str>   population type
        CP,f2,riself,bcsft
    --vcf <file> wgs_v4的vcf文件
    --p1 <str> 亲本1的id
    --p2 <str> 亲本2的id
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}
 mark=file(params.mark)
 trt =file(params.trt)
 map =file(params.map)
 anno_summary_file = file(params.anno)
 vcf = file(params.vcf)


process dataprepair {
    publishDir "${params.outdir}/01.dataprepair", pattern:"*"
    executor 'slurm'
    queue 'SANGERDEV'
    cpus 2
    memory "8G"
    input:
        file mark from mark
        file map from map
        file trt from trt
        file vcf from vcf
    output:
        file "total.rqtl.csv" into mark_file1,mark_file2
        file "total.map" into  map_file1,map_file2
        file "total.qtl.list" into qtl
        file "total.btl.list" into btl
        file "variant.diff.table" into variant_diff_table_ch
    script:
    if(params.popt == "CP"){
        """
		source ${params.env}
	ln -s ${map} total.map
        perl ${baseDir}/bin/map2rqtlCP.pl -l $mark -o total.rqtl.csv
        Rscript ${baseDir}/bin/trt.R -i $trt -o tmp -d total.rqtl.csv.order
		perl ${baseDir}/bin/trt.pl -i tmp.qtl.txt -o total.qtl --popt ${params.popt}
        perl ${baseDir}/bin/trt.pl -i tmp.btl.txt -o total.btl --popt ${params.popt}
        bcftools view -s ${params.p1},${params.p2} ${vcf} | bcftools query -H -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%TGT][\\t%AD]\\t%ANN\\n" | perl -ne 'next if(/#/);chomp;@a=split;@b=split(/\\,/,\$a[-1]);foreach \$b(@b){@c=split(/\\|/,\$b,-1);print join("\\t",@a[0..\$#a-1],@c),"\\n";}' > variant.table
        awk '{if((\$6 != "./." || \$5 != "./.") && \$5 != \$6) {safe++;print \$0}}' variant.table > variant.diff.table
        sed -i "1i CHROM\\tPOS\\tREF\\tALT\\t${params.p1}:GT\\t${params.p2}:GT\\t${params.p1}:AD\\t${params.p2}:AD\\tANN_Allele\\tANN_Annotation\\tAnnotation_impact\\tGene_Name\\tGene_ID\\tFeature_Type\\tFeature_ID\\tTranscript_BioType\\tRank\\tHGVS.c\\tHVGS.p\\tcDNA_stat\\tCDS_stat\\tAA_stat\\tDistance\\tERRORS" variant.diff.table
        """
    }else{
        """
		source ${params.env}
        perl ${baseDir}/bin/map2rqtl.pl -l $mark -m $map -o total.rqtl.csv
        ln -s ${map} total.map
        Rscript ${baseDir}/bin/trt.R -i $trt -o tmp
        perl ${baseDir}/bin/trt.pl -i tmp.qtl.txt -o total.qtl --popt ${params.popt}
        perl ${baseDir}/bin/trt.pl -i tmp.btl.txt -o total.btl --popt ${params.popt}
        bcftools view -s ${params.p1},${params.p2} ${vcf} | bcftools query -H -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%TGT][\\t%AD]\\t%ANN\\n" | perl -ne 'next if(/#/);chomp;@a=split;@b=split(/\\,/,\$a[-1]);foreach \$b(@b){@c=split(/\\|/,\$b,-1);print join("\\t",@a[0..\$#a-1],@c),"\\n";}' > variant.table
        awk '{if((\$6 != "./." || \$5 != "./.") && \$5 != \$6) {safe++;print \$0}}' variant.table > variant.diff.table
        sed -i "1i CHROM\\tPOS\\tREF\\tALT\\t${params.p1}:GT\\t${params.p2}:GT\\t${params.p1}:AD\\t${params.p2}:AD\\tANN_Allele\\tANN_Annotation\\tAnnotation_impact\\tGene_Name\\tGene_ID\\tFeature_Type\\tFeature_ID\\tTranscript_BioType\\tRank\\tHGVS.c\\tHVGS.p\\tcDNA_stat\\tCDS_stat\\tAA_stat\\tDistance\\tERRORS" variant.diff.table
        """
    }
} 


qtl.splitCsv(header:false,sep:'\t').groupTuple().set{qtl_list}
btl.splitCsv(header:false,sep:'\t').groupTuple().set{btl_list}
process qtl {
    publishDir "${params.outdir}/02.qtl", pattern:"*"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.pdf"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.result"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.png"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.pm.csv"
    executor 'slurm'
    queue 'SANGERDEV'
    cpus 2
    memory "8G"
    queueSize = 25
    input:
        file loc from mark_file1
        file map from  map_file1
        tuple qname,qtl_file from qtl_list
    output:
        file "*"
        tuple qname,"${qname}/${qname}.qtl-result.result"  into qtl_result
    script:
    if(params.popt == "CP"){
        """
		source ${params.env}
        Rscript ${baseDir}/bin/qtl-CP.R --map ${map} --loc ${loc} --trt ${qtl_file[0]} --out ${qname} --num 1000 --pvalue 0.01
        Rscript ${baseDir}/bin/manhattan.R -i ${qname}/${qname}.detail.result  -o ${qname}/${qname}
        """
    }else{
        """
		source ${params.env}
        Rscript ${baseDir}/bin/qtl-NOCP.R --map ${map} --loc ${loc} --trt ${qtl_file[0]} --out ${qname} --pop ${params.popt} --num 1000 --method cim --pvalue 0.01
        Rscript ${baseDir}/bin/manhattan.R -i ${qname}/${qname}.detail.result  -o  ${qname}/${qname}
        """
    }
} 
process btl {
    publishDir "${params.outdir}/03.btl", pattern:"*"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.pdf"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.result"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.png"
    //publishDir "${params.outdir}/07.result/${qname}/qtl", pattern:"*.pm.csv"
    executor 'slurm'
    queue 'SANGERDEV'
    cpus 2
    memory "8G"
    queueSize = 25
    cache 'lenient'
    input:
        file loc from mark_file2
        file map from  map_file2
        tuple qname,btl_file from btl_list
    output:
        file "*"
        tuple qname,"${qname}/${qname}.qtl-result.result"  into btl_result
    script:
    if(params.popt == "CP"){
        """
		source ${params.env}
        Rscript ${baseDir}/bin/qtl-CP.R --map ${map} --loc ${loc} --trt ${btl_file[0]} --out ${qname} --btl 1 --num 1000 --pvalue 0.001
        Rscript ${baseDir}/bin/manhattan.R -i ${qname}/${qname}.detail.result  -o ${qname}/${qname}
        """
    }else{
        """
		source ${params.env}
        Rscript ${baseDir}/bin/qtl-NOCP.R --map ${map} --loc ${loc} --trt ${btl_file[0]} --out ${qname}  --btl 1 --num 1000 --pvalue 0.001 --pop ${params.popt}
        Rscript ${baseDir}/bin/manhattan.R -i ${qname}/${qname}.detail.result  -o  ${qname}/${qname}
        """
    }
} 

regions=Channel.of().mix(qtl_result,btl_result)
process Orgdb {
    publishDir "${params.outdir}/00.orgdb", pattern:"*"
    executor 'slurm'
    queue "${params.queue}"
    cpus 2
    memory "8G"
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
    //publishDir "${params.outdir}/07.result/${method}/enrich", pattern:"${method}.region.gene.xls"
    //publishDir "${params.outdir}/07.result/${method}/enrich", pattern:"${method}.region.variant.xls"
    //publishDir "${params.outdir}/07.result/${method}/enrich/KEGG_result", pattern:"*KEGGenrichment*"
    //publishDir "${params.outdir}/07.result/${method}/enrich/GO_result", pattern:"*GOenrichment*"
    queue "${params.queue}"
    executor "slurm"
    cpus 2
    memory "8G"
    queueSize = 25
    input:
        tuple val(method),region from regions
        tuple path(term2gene),path(term2name) from term_ch
        file all_genes_list from all_gene_ch
        path LIB from EnrichDb
        file KEGG_anno from kegg_anno_ch
        file GO_anno from go_anno_ch
        file anno_summary from anno_summary_file
        file variant_diff_table from variant_diff_table_ch
    output:
        file "*"
    script:
    if(region.countLines() !=0){
    """
    #!/usr/bin/sh
        source ${params.env}
        mkdir ${method}
        cat ${region} | sed 1d | tr " " "\\t" | cut -f 6,7  > ${method}.region_list.txt
        Rscript ${baseDir}/bin/abstract_pop_genes.R --regionfile ${method}.region_list.txt --genefile ${all_genes_list} --outfile ${method}.genes_abstract.list
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$5}' > ${method}.degfile
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$6}' > ${method}.transcript.degfile
        Rscript ${baseDir}/bin/enrich.R --degfile ${method}.transcript.degfile --term2genefile ${term2gene} --term2namefile ${term2name} --outname ${method} --db ${LIB} --outdir ${method}/ --title ${method}
        Rscript ${baseDir}/bin/abstract_anno_summary.R --regionfile ${method}.region_list.txt --anno_summary_file ${anno_summary_file} --variant_table_file ${variant_diff_table} --outfile1 ${method}.region.gene.xls --outfile2 ${method}.region.variant.xls
        echo ${method}
    """
    }else{
    """
        #mkdir ${method}
        mkdir -p ${method}/KEGG_result
        touch ${method}/KEGG_result/noresult
        mkdir -p ${method}/GO_result
        touch ${method}/GO_result/noresult
        echo "no-region,please check!" > ${method}.error.log
    """
    }

}

def runcmd(cmd){
    println cmd
    proc = cmd.execute()
    def sout = new StringBuffer()
    proc.consumeProcessOutput(sout,sout)
    proc.waitFor()
    println sout.toString()
    assert !proc.exitValue()
}


workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    if(workflow.success){
        println "Workflow complete! ${workflow.duration}."
        //runcmd("cp -f ${projectDir}/bin/readme/qc结果说明文档.txt ${outdir}/qc/结果说明文档.txt")
        //runcmd("cp -f ${projectDir}/bin/readme/gs结果说明文档.txt ${outdir}/gs/结果说明文档.txt")
        runcmd("perl ${baseDir}/bin/gather.pl -i ${params.outdir} -o ${params.outdir}/")
        println "QTL complete!"
    } else {
        //println "Workflow failed! ${workflow.duration}."
    }
}
