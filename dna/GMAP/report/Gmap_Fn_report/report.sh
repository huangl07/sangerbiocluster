if [ -d "/mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report" ]; then

rm -rf /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report

fi

mkdir -p /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report

cp -r /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/src /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report

cp -r /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/css /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report

cat /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/title.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/workflow.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/1.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/2.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/3.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/4.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/5.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/6.*.rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/rmd/appendix.rmd > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report.rmd

sed -i -e 's/genome_chinese/西瓜/g' -e 's/genome_latin/Citrullus_lanatus/g' -e 's/sample_num/96/g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report.rmd 

/mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/bin/filter_rmd --rmd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report.rmd --format html --outfile /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_html.rmd

mkdir -p /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file

/mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/bin/get_project_info --project_contract MJ20190225074 --project_number FX2021083000054 --project_customer 蔡中全 --project_genome_chinese 西瓜 --project_genome Citrullus_lanatus --project_samplenum 96 --outfile /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/project_info.xls

sed -i 's/sample_Id/Mix_BuK/g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_html.rmd

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/project.info.xls   /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result.info.xls   /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file

cat /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/result.info.xls|while read line;do key=`echo "$line"|cut -f1`;value=`echo "$line"|cut -f2`;sed -i "s/${key}/${value}/g"  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_html.rmd;done 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/01.fastq_qc/rawdata_qc/qc.xls   /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/rawdata.xls

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/01.fastq_qc/cleandata_qc/atgc/Mix_BuK-1.atgc.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/01.fastq_qc/cleandata_qc/qual/Mix_BuK-1.qual.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/01.fastq_qc/cleandata_qc/qc.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

sed -i  's/#sampleID/Sample ID/g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/qc.xls && sed -i 's/#sampleID/Sample ID/g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/rawdata.xls &&sed -i '1d' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/*.atgc.xls

cd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file && for j in `ls *qual.xls`;do samples=`echo $j|sed 's/\.qual\.xls//g'`;less $j |awk '{print $1"	"$NF}'|sed '1d' > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/$samples.qual.sort.xls;Rscript /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/bin/ngsqc.r --base $samples.atgc.xls --qual $samples.qual.sort.xls --key $samples --od /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ ; done 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/03.map_stat/result.stat/Total.mapped.detail.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/03.map_stat/coverage/Mix_BuK.coverage.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/Total.mapped.detail.xls |cut -f1-4 > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/align_stat.xls 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/Total.mapped.detail.xls |awk -F "	" '{print$1"	"$8"	"$9"	"$6}' > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/coverage_sample.xls 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/Mix_BuK.coverage.xls |cut -f1|sort -u|grep -v 'sca'|sort -V > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/chrlist 

cd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file && for j in `ls *coverage.xls`;do samples=`echo $j|sed 's/\.coverage\.xls//g'`;Rscript /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/bin/genomeCoveragehorizontalArea.R --infile $j --idfile /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/chrlist --outfile /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/$samples.coverage --group.col 1 --x.col 2 --y.col 3 --x.lab Sequence-Position --y.lab AverageDepth-log2 --skip 0 --unit 100kb --log2 ; done

sed -i 's/^#//g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/align_stat.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/coverage_sample.xls

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/04.snp_indel/variant_stat/snp.stat  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/04.snp_indel/anno_stat/snp.stat  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/snp_anno.xls 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/snp.stat |cut -f1-7 > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/snp_stat.xls 

sed -i 's/^#//g' /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/snp_anno.xls 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/04.snp_indel/variant_stat/indel.stat  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/workflow_results/04.snp_indel/anno_stat/indel.stat  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/indel_anno.xls

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/indel.stat |cut -f1-7 > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/indel_stat.xls 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/01.genotype/pop.type.type.stat.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/01.genotype/pop.type.type.stat.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/03.gmap_evaluate/total.map.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/02.gmap/total.mapstat.xls|grep -v -E "：|:"|sed 's/#//g'  >/mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/total.mapstat.xls 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/03.gmap_evaluate/total.bin.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/03.gmap_evaluate/1.heatMap.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/03.gmap_evaluate/total.phy.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/03.gmap_evaluate/total.phy.spearman.xls|sed 's/^#//g'  >  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/total.phy.spearman.xls 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/04.qtl/yb.qtl.csv|sed 's/,/	/g' > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/t.qtl.xls 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/04.qtl/yb.scan.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/t.scan.png 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/05.anno/yb.gene.total|grep -E '^@|^#@'|sed -e 's/@//g' -e 's/#//g'|grep -E 'pos|chr4' > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/t.gene_anno.xls 

less /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/05.anno/yb.vcf.total|grep -E '^@|^#@'|sed  -e 's/@//g' -e 's/#//g'|grep -E 'pos|chr4' > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/region.anno.xls 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/06.enrich/GO_result/genes_ORA_GOenrichment_bubble.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/06.enrich/GO_result/genes_GOenrichment.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/06.enrich/KEGG_result/genes_ORA_KEGGenrichment_bubble.png  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cp /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Result/result/06.enrich/KEGG_result/genes_ORA_KEGGenrichment.xls  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file/ 

cd /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/file 

rm *pdf

Rscript  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/Gmap_Fn_report/bin/rmarkdown --rmd  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_html.rmd --format html --outfile /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_raw.html

grep -v -E "Warning in stringi|integer vector|Warning in instance|^##"  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_raw.html > /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/FX2021083000054_report.html

if [ -s /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/FX2021083000054_report.html ];then

    rm  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report.rmd

    rm  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_raw.html

    rm  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/report_html.rmd

    rm -rf  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/css

    rm -rf  /mnt/ilustre/users/jiawen.ma/development/Result_report/Gmap_Fn_report/FX2021083000054_report/src

else

exit 1

fi
