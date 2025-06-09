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
params.deepbsa=false
params.noparent=false
params.env="~/app/bioinfo/dna/new.rc_1"
params.queue="SANGERDEV"
params.qtgseq=false
params.boot=false
params.minIndex=0.3
params.minmarker=10
params.ems=false
def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    --vcf   <file>  input vcf file [forced]
    --outdir   <dir>   output dir [forced]
    --chr   <file>  chr list for draw [forced]
    ##############################################################
    --gff   <file>  gff file [forced]
    --anno_summary    <file>  anno_summary_file [forced] eg tmp/02.reference/anno.summary
    ##############################################################
    --group <file>   group file  [forced]
        mbid must be given

        wpid    flag    phdep   pldep
        mpid    flag    phdep   pldep

    --popt  <str>   population [forced]
    ###############################################################
    --bulksize  <num>   bulk size default 30
    --winsize <num>     windows size default 1000000
    --stepsize <num>    step size default 1000
    --minmarker <num>   report region which include at least this number of markers default 10
    ###############################################################
    --boot <bool>   use bootstrap method default false
    --bootstrap <num>   bootstrap number default 1000
    --pvalue    <num>   pvalue
        ridit for pvalue [ 0.0001]
    ################################################################
    --grade <bool>  do grade for multi bulk than 2
    --noparent <bool> do no parents or abs index
    --mutmap <bool> do mutmap
    --minIndex <float> mutmap min index threshold
    --ems <bool> ems filter
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
        if (params.ems && params.mutmap)
            """
            source ${params.env}
            perl ${baseDir}/bin/EMS_select.pl --vcf ${vcf} --out ems.table --group ${group}
            cut -f 1-3 ems.table > ems.list
            bcftools view --threads 2 -T ems.list -Oz -o ems.vcf.gz ${vcf}
            perl ${baseDir}/bin/vcf2table.pl --vcf ems.vcf.gz --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
            perl ${baseDir}/bin/bsa_calc.pl --table pop.table --group ${group} --popt ${params.popt} -out pop.index --minindex ${params.minIndex}
            perl ${baseDir}/bin/anno.pl --table pop.anno --anno ${anno} -out pop.final.anno
            perl ${baseDir}/bin/bsa_stat_filtered_pop_table.pl --infile pop.table --outfile snp_indel_gene.stat --chrfile ${chr} --anno ${anno}
            """
        else
            """
            source ${params.env}
            perl ${baseDir}/bin/vcf2table.pl --vcf ${vcf} --out pop.table --group ${group} -popt ${params.popt} -vtype ALL
            perl ${baseDir}/bin/bsa_calc.pl --table pop.table --group ${group} --popt ${params.popt} -out pop.index --minindex ${params.minIndex}
            perl ${baseDir}/bin/anno.pl --table pop.anno --anno ${anno} -out pop.final.anno
            perl ${baseDir}/bin/bsa_stat_filtered_pop_table.pl --infile pop.table --outfile snp_indel_gene.stat --chrfile ${chr} --anno ${anno}
            """
    }

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
            sed 1d ${pop_table} | awk  -v OFS='\\t' '{print \$1,\$2,\$2}' > pos_list.txt
            sed -i "1d" pos_list.txt
            tabix ${vcf}|| tabix -C ${vcf}
            python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list --pop_info pop.list
            bcftools view -S mix.list -R pos_list.txt ${vcf} -o filter.vcf -Ov --threads 8
            bcftools view -S pop.list -R pos_list.txt ${vcf} -o filter.vcf.gz -Oz --threads 8
        """
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
        def abs = params.noparent ? "--abs": ""
        def mutmap = params.mutmap ? "--mutmap": ""
        def method = (params.winsize <= 1000 && params.stepsize <= 1000) ? "num": "bp"
        if(params.boot){
            """
            source ${params.env}
            Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize} --method ${method} --minmarker ${params.minmarker} ${abs}
            Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue} ${mutmap}
            Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --qcol CI --xlab chromosome --ylab "delta-index" ${mutmap} ${abs} --method ${method}
            Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number ${params.minmarker}
            echo index_winsize:${params.winsize} > index.params.log
            echo index_stepsize:${params.stepsize} >> index.params.log
            echo index_bootstrap:${params.bootstrap} >> index.params.log
            echo index_p:${params.pvalue} >> index.params.log
            echo minmarker:${params.minmarker} >> index.params.log
            """
        } else {
            """
            source ${params.env}
            Rscript ${baseDir}/bin/slidingwin-index_test.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize} --method ${method} --minmarker ${params.minmarker} ${abs}
            Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.sliding.detail --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue} ${mutmap} --noboot
            Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result pop.bootstrap.result --chr ${chr} --output pop.index --pcol delta --lcol slidingD --threshold ${params.quantile} --xlab chromosome --ylab "delta-index" ${mutmap} ${abs} --method ${method}
            Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --CI CI --outfile pop.index.region --number ${params.minmarker}
            echo index_winsize:${params.winsize} > index.params.log
            echo index_stepsize:${params.stepsize} >> index.params.log
            echo index_p:${params.pvalue} >> index.params.log
            echo index_bootstrap:0 >> index.params.log
            echo index_minmarker:${params.minmarker} >> index.params.log
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
        def mutmap = params.mutmap ? "--mutmap": ""
        def abs = params.noparent ? "--abs": ""
        """
            source ${params.env}
            Rscript ${baseDir}/bin/loess-index-snow.R --infile ${table_file} --out pop.index.loess ${mutmap}
            Rscript ${baseDir}/bin/bootstrap-index.R --infile pop.index.loess --bulksize ${params.bulksize} --outfile pop --bootstrap ${params.bootstrap} --popstruc ${params.popt} --qvalue ${params.pvalue} ${mutmap} --noboot ${abs}
            Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result pop.bootstrap.result --pcol delta --lcol slidingD --xlab chromosome --ylab "loess" --threshold ${params.quantile}  --output loess --chr ${chr} ${mutmap} ${abs}
            Rscript ${baseDir}/bin/region.R --infile pop.bootstrap.result --ccol X.chr --pos pos --loess slidingD --quantile ${params.quantile} --outfile pop.loess.region --number ${params.minmarker}
            echo loess_q:${params.quantile} > loess.params.log
            echo loess_minmarker:${params.minmarker} >> loess.params.log
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
            def method = (params.winsize <= 1000 && params.stepsize <= 1000) ? "num": "bp"
            """
                source ${params.env}
                Rscript ${baseDir}/bin/G-calc.R --infile ${index_file} --outfile pop.Gprime --winsize ${params.winsize} --method ${method} --threshold ${params.pvalue}
                Rscript ${baseDir}/bin/manhattan_bsa_v2_index_quantile.R --result pop.Gprime --pcol G --lcol Gprime --threshold ${params.quantile} --xlab chromosome --ylab "Gprime" --output Gprime  --chr ${chr}
                Rscript ${baseDir}/bin/region.R --infile pop.Gprime --ccol X.chr --pos pos --loess Gprime --quantile ${params.quantile} --outfile pop.Gprime.region --number ${params.minmarker}
                echo gprime_winsize:${params.winsize} > gprime.params.log
                echo gprime_q:${params.quantile} >> gprime.params.log
                echo gprime_minmarker:${params.minmarker} >> gprime.params.log
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
        def method = (params.winsize <= 1000 && params.stepsize <= 1000) ? "num": "bp"
        """
            source ${params.env}
            Rscript ${baseDir}/bin/slidingwin-index.R --infile ${index_file} --outfile pop --winsize ${params.winsize} --stepsize  ${params.stepsize}  --method ${method}
            Rscript ${baseDir}/bin/manhattan_bsa_v2.R --result pop.sliding.result --chr ${chr} --output pop.ED  --pcol slidingED --threshold ${params.quantile} --xlab chromosome --ylab "ED"
            Rscript ${baseDir}/bin/region.R --infile pop.sliding.detail --ccol X.chr --pos pos --loess slidingED --quantile ${params.quantile} --outfile pop.ED.region --number ${params.minmarker}
            uniq pop.sliding.result
            echo ed_winsize:${params.winsize} > ed.params.log
            echo ed_stepsize:${params.stepsize} >> ed.params.log
            echo ed_q:${params.quantile} >> ed.params.log
            echo ed_minmarker:${params.minmarker} >> ed.params.log
        """
        }
        if(params.qtgseq){
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
        }
        if(params.deepbsa){
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
            Rscript ${baseDir}/bin/region.R --infile DL_values.txt --ccol X.chr --pos pos --loess loess --quantile ${params.quantile} --outfile pop.deepBSA_DL.region --number ${params.minmarker}
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
            Rscript ${baseDir}/bin/region.R --infile K_values.txt --ccol X.chr --pos pos --loess loess --quantile ${params.quantile} --outfile pop.deepBSA_K.region --number ${params.minmarker}
            """
            }
        }
    }


    regions = Channel.of().mix(loess_region,index_region)
    if(!params.mutmap){
        regions = regions.mix(ED_region,Gprime_region)
    }
    if(params.qtgseq){
        regions = regions.mix(QTG_region)
    }
    if(params.deepbsa){
        regions = regions.mix(deepBSA_DL_region,deepBSA_K_region)
    }
    regions.multiMap{it -> r1: r2: it}.set {regions_ch}

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
            Rscript ${baseDir}/bin/region-ridit.R --infile pop.denoise.result --ccol X.chr --pos pos --loess logP --CI CI --outfile pop.ridit.region --number ${params.minmarker} --winsize ${params.winsize}
        """
    }

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
                tabix ${vcf} || tabix -C ${vcf}
                python3 ${baseDir}/bin/abstract_mix_info.py --group_file ${group} --mix_info mix.list
                bcftools view -S mix.list -R pos_list.txt ${vcf} -o deepBSA.vcf -Ov --threads 8
            """
    }
}

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
    queue "${params.queue}"
    executor "slurm"
    tag "${method}"
    cpus 4
    memory "32G"
    input:
        tuple val(method),region from regions_ch.r1
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
        tabix ${vcf_file} || tabix -C ${vcf_file}
        mkdir ${method}
        Rscript ${baseDir}/bin/abstract_pop_genes.R --regionfile ${region} --genefile ${all_genes_list} --outfile ${method}.genes_abstract.list
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$5}' > ${method}.degfile.tmp
        cat ${method}.genes_abstract.list | awk 'NR>1{print \$6}' > ${method}.transcript.degfile.tmp
        cut -f5,6 ${method}.genes_abstract.list  > geneID_transcriptID
        csvtk cut -f GeneName,GeneID -d '\t' -T pop.final.anno >geneNAME_geneID
        csvtk merge -L --na -- -Tt -f "gene_id;GeneID" geneID_transcriptID geneNAME_geneID >geneID_transcriptID_geneNAME 
        Rscript ${baseDir}/bin/enrich.R --degfile ${method}.transcript.degfile.tmp --term2genefile ${term2gene} --term2namefile ${term2name} --outname ${method} --db ${LIB} --outdir ${method}/ --title ${method} --genename geneID_transcriptID_geneNAME
        echo -e "Gene_id\\tTranscript_id\\tChr\\tPos1\\tPos2\\tNR_ID\\tNR_ANNO\\tUni_ID\\tUni_ANNO\\tKEGG_ID\\tKEGG_ANNO\\tGO_ID\\tGO_ANNO\\tEggNOG_ID\\tEggNOG_ANNO\\tPfamAccession\\tPfamAnno" > ${method}.degfile
        cp ${method}.degfile ${method}.transcript.degfile
        grep -wf ${method}.degfile.tmp ${anno_summary} | sed "s/:/\\t/g" >> ${method}.degfile
        grep -wf ${method}.transcript.degfile.tmp ${anno_summary} | sed "s/:/\\t/g" >> ${method}.transcript.degfile
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

process extract_region_from_vcf{
    publishDir "${params.outdir}/07.extract_vcf", pattern: "*.xls"
    queue "${params.queue}"
    executor "slurm"
    tag "${method}"
    cpus 4
    memory "32G"
    input:
        tuple val(method),region from regions_ch.r2
        file vcf_file from pop_filtered_vcf_ch3
        file anno_summary_file from anno_summary_file
    output:
        tuple val(method), file("${method}.extract.vcf") into extract_vcf_ch, extract_vcf_ch2
        tuple val(method), file("${method}.extract.vcf.gz") into extract_gz_vcf_ch
        file("${method}.anno.xls")
    script:
        """
        source ${params.env}
        tabix -f ${vcf_file} || tabix -C ${vcf_file}
        bcftools view --threads 4 -R ${region} ${vcf_file} -o ${method}.extract.vcf.gz
        gunzip -k ${method}.extract.vcf.gz
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%TYPE\\t[%TGT\\t][%AD\\t]%ANN\\n' -H ${method}.extract.vcf.gz | perl ${baseDir}/bin/vcf2table_wgs.pl > ${method}.table
        Rscript ${baseDir}/bin/anno_summary_to_variant_table.R -a ${anno_summary_file} -v ${method}.table -o ${method}.anno.xls
        """
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
