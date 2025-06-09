process params.stat{
    publishDir "${params.outdir}/05.enrich", pattern:"*"
    queue "DNA"
    executor "slurm"
    input:
    output:
        file "*"
    script:
    if(region.countLines() != 0){
    """
        mkdir ${method}
        cd ${method}
        Rscript ${baseDir}/bin/abstract_pop_genes.R --regionfile ${region} --genefile ${all_genes_list} --outfile ../genes_abstract.list
        cat genes_abstract.list | awk 'NR>1{print \$4}' > degfile
        Rscript ${baseDir}/bin/enrich.R --degfile degfile --term2genefile ../${term2gene} --term2namefile ../${term2name} --outname ${method} --db ../${LIB}
        perl ${baseDir}/bin/extract_region.gene.eff.pl -gff ../${gff} -region ${region} -table ../${table} -out ${method}
        Rscript ${baseDir}/bin/merge_annotation.R --regionfile ${region} --genefile genes_abstract.list --gofile ../${GO_anno} --keggfile ../${KEGG_anno} --outfile region_gene.txt
    """
    }
}