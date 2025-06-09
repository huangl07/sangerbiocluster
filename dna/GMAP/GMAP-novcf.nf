#!/usr/bin/env nextflow
params.out = "demo"
params.popt="F2"
params.segment=0.5
params.missing=0.3
params.Pdep=10
params.Odep=2
params.nchro=1
params.help = false
params.chr_only=false
params.vcf = "published/data/04.snpIndel/finalvcf/pop.snpindel.final.vcf.gz"
def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    --filterd <file> marker file
    --popt  <str>   population type CP/BCi/Fi/RILx/
    --out   <dir>   output dir
    --p1    <str>   input p1 IDs
    --p2    <str>   input p2 IDs
    --Pdep  <num>   input parents depth 10
    --Odep  <num>   input offspring depth 5
    --segment   <num>   segment 0.05 more bigger is less accuracy
    --SSR   <file>  input SSR merge file
    --missing   <num>   missing 0.5
    --chr   <file>  chr.list
    --chr_only

    --nchro <num>   chr number
     """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

filterd=file(params.filterd)
chr=file(params.chr)
process markersplit{
    publishDir "${params.out}/02.markersplit", pattern:"*"
    cpus 8
    memory "32G"
    input:
        file marker from filterd
        file chrfile from chr
    output:
        file "scaf.list" into scaflist
        file "split/compare.sh" into compare_worksh,compare
        path "split"
        file "*"
    script:

    if(params.popt == "CP"){
      if(params.chr){
            """
            source ~/app/bioinfo/dna/new.rc
            perl ${baseDir}/bin/markersplitCP.pl -input ${marker} -fOut scaf.list -dOut split -chr ${chrfile}
            """
        }else{
            """
            source ~/app/bioinfo/dna/new.rc
            perl ${baseDir}/bin/markersplitCP.pl -input ${marker} -fOut scaf.list -dOut split
            """
        }
    }else{
        if(params.chr){
            """
            source ~/app/bioinfo/dna/new.rc
            perl ${baseDir}/bin/markersplitNOCP.pl -input ${marker} -fOut scaf.list -dOut split -chr ${chrfile} 
            """
        }else{
            """
            source ~/app/bioinfo/dna/new.rc
            perl ${baseDir}/bin/markersplitNOCP.pl -input ${marker} -fOut scaf.list -dOut split 
            """
        }
    }
}


scaflist.splitCsv(header:false,sep:'\t').groupTuple().set{scaf_list}

process binmap{
    publishDir "${params.out}/03.bin", pattern:"*"
    cpus 8
    memory "32G"
    input:
        tuple val(sca),val(len),scaf from scaf_list
    output:
        file "*.bin.marker" into binfile1,binfile2
        file "*"
    script:
        if(params.popt == "CP"){
            """
            source ~/app/bioinfo/dna/new.rc
            ln -s ${scaf[0]} ${sca}.bin.marker
            """
        }else{
            if(file("${scaf[0]}").countLines() > 500){
            """
            source ~/app/bioinfo/dna/new.rc            
            less -S ${scaf[0]} |sort -n -k 2 > ${sca}.sort.marker
            snpbinner crosspoints -i ${sca}.sort.marker -o ${sca}.cross -r 0.001
            snpbinner bins -i ${sca}.cross -o ${sca}.bins -l 1000
            perl ${baseDir}/bin/convert2MSTmap.pl -input ${sca}.bins -output ${sca}.bin.marker -popt ${params.popt} --chr ${sca} --marker ${scaf[0]}
            """
            }else{
            """
            source ~/app/bioinfo/dna/new.rc
            less -S ${scaf[0]}|sed 's/\t/-/'|sed 's/chr-pos/MarkerID/g' > ${sca}.tmp.marker
            less -s ${sca}.tmp.marker|perl -ne 'chomp;if(/MarkerID/){print \$_,"\n";}else{(\$a,\$b)=split(/\t/,\$_,2);\$b=~s/-/U/g;\$b=~s/h/X/g;print \$a,"\t",\$b,"\n";}' > ${sca}.bin.marker
            """
            }
        }
}

onlychr=0
if(compare.splitText(by:1).toList().size().val < 1 || params.chr_only){
   println "haha"
   onlychr=1
   process group_by_chr{
        cpus 2
        memory "32G"
        publishDir "${params.out}/06.premapping", pattern:"*"
        input:
            file bins from binfile1.collect()
        output:
            file "*.lg" into grouping_file
            file "*"
        script:
        """
        source ~/app/bioinfo/dna/new.rc
        cat *.marker|grep -f ${chr} > total.markers
        perl ${baseDir}/bin/linkage_by_ref.pl  -i total.markers -o  total.lg
        """
    }

}else{
    compare_worksh.splitText(by:1).set{compare_para}
    println "test"
    process mlodcalc{
        publishDir "${params.out}/04.calc", pattern:"*"
        cpus 2
        memory "32G"
        input:
            val para from compare_para
            file binfile from binfile1.collect()
        output:
            file "*.mlod" into mlod_file
            file "*"
        script:
        """
        perl  ${baseDir}/bin/calculateMLOD.pl  -popt ${params.popt} ${para}
        """
    }
    process grouping{
        publishDir "${params.out}/05.group", pattern:"*"
        cpus 2
        memory "32G"
        input:
            file mlod from mlod_file.collect()
            file chrfile from chr
        output:
            file "*.lg" into grouping_file
            file "*"
        script:
            if(params.chr){
                """
                            source ~/app/bioinfo/dna/new.rc

                cat *.mlod > Total.mLOD
                python3  ${baseDir}/bin/count_mlod.py -input Total.mLOD -output linkage.mlod.csv
                perl ${baseDir}/bin/linkage_by_ref-mlod.pl  -i linkage.mlod.csv -o lg.lg -t 5 -c ${chrfile}
                """
            }else{
                """
                            source ~/app/bioinfo/dna/new.rc

                cat *.mlod > Total.mLOD
                python3  ${baseDir}/bin/count_mlod.py -input Total.mLOD -output linkage.mlod.csv
                perl ${baseDir}/bin/linkage_by_mlod.pl  -i linkage.mlod.csv -k lg -d ./ -n ${params.nchro} -minGroup 1 -b 3 -e 20
                """
            }
    }
}

    process premapping{
        publishDir "${params.out}/06.premapping", pattern:"*"
        cpus 2
        memory "32G"
        input:
            file lg from grouping_file
            file bins from binfile2.collect()
        output:
            file "linkagegroups/lg.list" into lglist
            file "linkagegroups/*.marker" into lgmarkers1,lgmakrers2
            file "head" into head
            file "*"
        script:
        if(params.popt == "CP"){
            """
                cat *.marker > total.markers
                perl ${baseDir}/bin/splitbyLG-CP.pl  -l ${lg} -i total.markers -d linkagegroups/
                head -n 1 total.markers|sed 's/type/type\tphase/' > head
            """
        }else{
            if(onlychr > 0){
            """
            source ~/app/bioinfo/dna/new.rc
            cat *.marker > total.markers
            perl ${baseDir}/bin/splitbyLG-NOCP.pl  -l ${lg} -i total.markers -d linkagegroups/ --chr 1
            head -n 1 total.markers > head

            """
            }else{
            """
            source ~/app/bioinfo/dna/new.rc
            cat *.marker|grep MarkerID|uniq > head
            cat *.marker|sort -r|less -S |uniq|grep -v pos > total.markers.tmp
            cat head total.markers.tmp > total.markers
            perl ${baseDir}/bin/splitbyLG-NOCP.pl  -l ${lg} -i total.markers -d linkagegroups/
            """
            }
        }
    }


lglist.splitCsv(header:false,sep:'\t').groupTuple().set{lg_list}

process mapping{
    publishDir "${params.out}/07.mapping", pattern:"*"
    publishDir "${params.out}/09.result", pattern:"*.pdf"
	publishDir "${params.out}/09.result",pattern:"*.png"
    cpus 2
    memory "32G"
    input:
        tuple lg,marker from lg_list
    output:
        file "${lg}.result.map" into mapfile
        file "${lg}.result.csv" into csv_result
        file "*"
    script:
    if(params.popt == "CP"){
        if(params.chr_only){
            """
            source ~/app/bioinfo/dna/new.rc
            crosslink_group --inp=${marker[0]} --outbase=${lg}. --knn=30 --matpat_lod=20 --matpat_weights=01P07 --min_lod=4 --redundancy_lod=20 --redun=${lg}.redun.list
            ln -s  `wc -l ${lg}.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.pri.draw.loc
            mcl  ${lg}.redun.list --abc  -I 2.0
            less -S  out.{$lg}.redun.list.I20|perl -ne 'chomp;@a=split;print join("\\n",@a[1..\$#a]),"\\n"' > ${lg}.redun.list; #这里可以考虑优化
            grep -wvf ${lg}.redun.list ${lg}.pri.draw.loc > ${lg}.draw.loc
            cut -f 1 -d " " ${lg}.draw.loc|perl -ne 'chomp;@a=split(/-/,\$_);print \$_,"\\t",\$a[1],"\\n";' > ref.map
            perl ${baseDir}/bin/smooth-CP.pl -l ${lg}.draw.loc -k ${lg} -d ./ -m ref.map   -ami 0.8
            crosslink_group --inp=${lg}.correct.loc --outbase=${lg}.correct. --knn=45 --map=${lg}.correct.
            ln -s  `wc -l ${lg}.correct.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.csv
            ln -s  `wc -l ${lg}.correct.*.map|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.tmp.map
            less -S ${lg}.tmp.map|grep -v group|perl -ne 'if(/group/){print "group\\t${lg}\\n"}else{@a=split;print join("\\t",@a),"\\n"}' > ${lg}.tmp.format.map
            Rscript ${baseDir}/bin/genericMap.R --input ${lg}.tmp.format.map --output ${lg}.result.map
            less -S ${lg}.result.map|perl -ne 'chomp;if(/group/){print \$_,"\n"}else{@a=split;print join("\\t",\$a[0],\$a[3]),"\n"}' > ${lg}.map
            sed -i '1 i group\\t${lg}' ${lg}.map
            perl ${baseDir}/bin/marker2joinmap.pl -i ${lg}.result.csv -o ${lg}
            Rscript ${baseDir}/bin/drawmap.R --mark ${lg}  --out ./ --pop ${params.popt}
            """
        }else if(params.chr){
            """
            source ~/app/bioinfo/dna/new.rc
            crosslink_group --inp=${marker[0]} --outbase=${lg}. --knn=30 --matpat_lod=20 --matpat_weights=01P07 --min_lod=4 --redundancy_lod=20 --redun=${lg}.redun.list
            ln -s  `wc -l ${lg}.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.pri.draw.loc
            #cut -f 1 -d " " ${lg}.draw.loc|perl -ne 'chomp;@a=split(/-/,\$_);print \$_,"\\t",\$a[1],"\\n";' > ref.map
            mcl  ${lg}.redun.list --abc  -I 2.0
            less -S  out.${lg}.redun.list.I20|perl -ne 'chomp;@a=split;print join("\\n",@a[1..\$#a]),"\\n"' > ${lg}.redun.list; # 这里可以考虑优化
            grep -wvf ${lg}.redun.list ${lg}.pri.draw.loc > ${lg}.draw.loc #
            cut -f 1 -d " " ${lg}.draw.loc|perl -ne 'chomp;@a=split(/-/,\$_);print \$_,"\\t",\$a[1],"\\n";' > ref.map
			perl ${baseDir}/bin/smooth-CP.pl -l ${lg}.draw.loc -k ${lg} -d ./ -m ref.map   -ami 0.8
            crosslink_group --inp=${lg}.correct.loc --outbase=${lg}.correct. --knn=45 --map=${lg}.correct.
            ln -s  `wc -l ${lg}.correct.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.csv
            ln -s  `wc -l ${lg}.correct.*.map|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.tmp.map
            less -S ${lg}.tmp.map|grep -v group|perl -ne 'if(/group/){print "group\\t${lg}\\n"}else{@a=split;print join("\\t",@a),"\\n"}' > ${lg}.tmp.format.map
            Rscript ${baseDir}/bin/genericMap.R --input ${lg}.tmp.format.map --output ${lg}.result.map
            less -S ${lg}.result.map|perl -ne 'chomp;if(/group/){print \$_,"\n"}else{@a=split;print join("\\t",\$a[0],\$a[3]),"\n"}' > ${lg}.map
            sed -i '1 i group\\t${lg}' ${lg}.map
            perl ${baseDir}/bin/marker2joinmap.pl -i ${lg}.result.csv -o ${lg}
            Rscript ${baseDir}/bin/drawmap.R --mark ${lg}  --out ./ --pop ${params.popt}
            """
        }else{
            """
            source ~/app/bioinfo/dna/new.rc
            crosslink_group --inp=${marker[0]} --outbase=${lg}. --knn=30 --matpat_lod=20 --matpat_weights=01P07 --min_lod=4 --redundancy_lod=20
            ln -s `wc -l ${lg}.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.draw.loc
            crosslink_map --inp=${lg}.draw.loc --map=${lg}.primary.map --out=${lg}.result.csv --ga_gibbs_cycles=10 --ga_iters=300000  --ga_max_hop=1.0 --ga_max_mvseg=1.0 --ga_max_mvdist=1.0 --ga_max_seg=1.0 --gibbs_samples=500  --gibbs_burnin=20 --gibbs_min_prob_1=0.1 --gibbs_min_prob_2=1
            perl ${baseDir}/bin/smooth-CP.pl -l ${lg}.result.csv -k ${lg} -d ./ -m ${lg}.primary.map -win 30 -ami 0.8
            crosslink_group --inp=${lg}.correct.loc --outbase=${lg}.correct.
            crosslink_map --inp=${lg}.correct.000.loc --map=${lg}.correct.000.map --ga_gibbs_cycles=10 --ga_iters=300000  --ga_max_hop=1.0 --ga_max_mvseg=1.0 --ga_max_mvdist=1.0 --ga_max_seg=1.0 --gibbs_samples=500  --gibbs_burnin=20 --gibbs_min_prob_1=0.1 --gibbs_min_prob_2=1
            ln -s  `wc -l ${lg}.correct.*.loc|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.result.csv
            ln -s  `wc -l ${lg}.correct.*.map|sort -nr|head -n2|tail -n 1|awk '{print \$2}'` ${lg}.tmp.map
            less -S ${lg}.tmp.map|grep -v group|perl -ne 'chomp;if(/group/){print "group\\t${lg}\\n"}else{@a=split;print join("\\t",@a),"\\n"}' > ${lg}.tmp.format.map
            Rscript ${baseDir}/bin/genericMap.R --input ${lg}.tmp.format.map --output ${lg}.result.map
            less -S ${lg}.result.map|perl -ne 'chomp;if(/group/){print \$_,"\n"}else{@a=split;print join("\\t",\$a[0],\$a[3]),"\n"}' > ${lg}.map
            sed -i '1 i group\\t${lg}' ${lg}.map
            perl ${baseDir}/bin/marker2joinmap.pl -i ${lg}.result.csv -o ${lg}
            Rscript ${baseDir}/bin/drawmap.R --mark ${lg}  --out ./ --pop ${params.popt}
            """
        }
    }else{
        """
        source ~/app/bioinfo/dna/new.rc
        Rscript ${baseDir}/bin/Asmap.R  --binfile ${marker[0]} --output ${lg}.result.map --popt ${params.popt} --lg ${lg}
        """
    }
}

println(mapfile)
process mapEvaluate{
    publishDir "${params.out}/08.evaluate", pattern:"*"
    publishDir "${params.out}/09.result", pattern:"total.result.csvs"
    publishDir "${params.out}/09.result", pattern:"total*map*"
    publishDir "${params.out}/09.result", pattern:"total.filtered.stat"
    publishDir "${params.out}/09.result", pattern:"total.mapstat"
    publishDir "${params.out}/09.result", pattern:"total.phy.spearman.xls"
    publishDir "${params.out}/09.result", pattern:"*.pdf"
    publishDir "${params.out}/09.result", pattern:"*.png"
    input:
        file map from mapfile.collect()
        file lg from lgmakrers2.collect()
        file csv from csv_result.collect()
        file head from head
    output:
        file "*"
    script:
    if(params.popt == "CP"){
        """
            source ~/app/bioinfo/dna/new.rc
            ln -s  ${stat} total.filtered.stat
            cat ${head} *.result.csv > total.result.csvs
            perl ${baseDir}/bin/marker2joinmap.pl -i total.result.csvs -o total
            perl ${baseDir}/bin/map-gather.pl -i ./ -o ./
            perl ${baseDir}/bin/marker2phase.pl -i total.result.csvs -o total.result.phase -m total.sexAver.map
            perl ${baseDir}/bin/mapEstimate.pl -i total.sexAver.map -o total.mapstat
            perl ${baseDir}/bin/mapEstimate.pl -i total.male.map -o total.male.mapstat
            perl ${baseDir}/bin/mapEstimate.pl -i total.female.map -o total.female.mapstat
            perl ${baseDir}/bin/drawAligmentRalationMap.pl -m total.sexAver.map -o ./ -k total.phy
            perl ${baseDir}/bin/markerinfo.pl -map total.sexAver.map -input total.result.csvs --pop ${params.popt} --out total
            less -S total.marker.info |cut -f 1,2,3|sed 's/\\t/,/g' > total.map.draw
            Rscript ${baseDir}/bin/plotmaps.R --mark total.map.draw --out total.map
	        Rscript ${baseDir}/bin/drawbinCP-sexAver.R --mark total.result.phase.draw  --out total.sexAver.bin
            #Rscript ${baseDir}/bin/drawmap.R --mark total  --out ./ --pop ${params.popt}
        """
    }else{
        """
        source ~/app/bioinfo/dna/new.rc
         cat  *.map > total.maps
         cat *.marker|grep -v MarkerID > tmp.markers
         cat ${head} tmp.markers > total.markers
          less -S *.result.csv|grep Genotype|uniq > tmp.head
         cat *.csv|sort|uniq|grep -v Genotype > tmp.csvs
         cat tmp.head tmp.csvs > total.csvs
         rm tmp.*
         perl ${baseDir}/bin/mapEstimate.pl -i total.maps -o total.mapstat
        perl ${baseDir}/bin/markerinfo.pl -map total.maps -input total.markers --pop ${params.popt} --out total
        Rscript ${baseDir}/bin/drawmap.R --mark total.csvs  --out ./ --pop ${params.popt}
        Rscript ${baseDir}/bin/plotmaps.R --mark total.csvs --out total.map
        Rscript ${baseDir}/bin/drawbinNOCP.R --mark total.csvs  --out ./ --name total.bin
        perl ${baseDir}/bin/drawAligmentRalationMap.pl -m total.maps -o ./ -k total.phy
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

