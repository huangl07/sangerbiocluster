#!bin/bash

workdir=`pwd`
data_realease_dir="$workdir/data_release"
result_dir="$workdir/result/"

cd $workdir
if [ ! -d $data_realease_dir ];then
    mkdir -p $data_realease_dir
fi

cd $data_realease_dir
mkdir 01.vcf2table 02.index 03.loess  04.enrich

# 01.vcf2table结果整理
cut -f 3 $result_dir/01.vcf2table/pop.index | sed 's/delta/index/'> tmp2.txt
cat $result_dir/01.vcf2table/pop.table |cut -f 1-5 >tmp1.txt
awk -F "\t" '{for (i=6;i<=NF;i++)printf("%s\t", $i);print ""}' $result_dir/01.vcf2table/pop.table >tmp3.txt
paste tmp1.txt tmp2.txt tmp3.txt > ./01.vcf2table/pop.table
rm tmp1.txt tmp2.txt tmp3.txt
cp $result_dir/01.vcf2table/snp_indel_gene.stat ./01.vcf2table/snp_indel_gene.stat.xls
cp $result_dir/01.vcf2table/pop.final.anno ./01.vcf2table/pop.final.anno
awk -v OFS='\t' '{$NF="";$(NF-1)="";$(NF-2)="";print $0}' ./01.vcf2table/pop.final.anno > ./01.vcf2table/pop.final.anno.xls
rm -rf ./01.vcf2table/pop.final.anno

# 02.index-slid结果整理
awk '{print $1"\t"$4"\t"$5"\t"$6}' $result_dir/02.index-slid/pop.bootstrap.result > ./02.index/pop.index.bootstrap.result.xls
cp $result_dir/02.index-slid/pop.index*.png ./02.index
cp $result_dir/02.index-slid/pop.index*.pdf ./02.index
if [ -d $result_dir/06.enrich/index ];then
  cp $result_dir/06.enrich/index/index.all.table ./02.index/index.all.table
  awk -v OFS='\t' '{$NF="";$(NF-1)="";$(NF-2)="";print $0}' ./02.index/index.all.table > ./02.index/index.all.table.xls
  rm -rf ./02.index/index.all.table
  cp $result_dir/06.enrich/index/index.region_stat.xls ./02.index/index.region_stat.xls
fi

# 03.loess结果整理
awk '{print $1"\t"$4"\t"$5"\t"$6}' $result_dir/05.loess/pop.bootstrap.result >./03.loess/pop.loess.bootstrap.result.xls
cp $result_dir/05.loess/*.png ./03.loess
cp $result_dir/05.loess/*.pdf ./03.loess
if [ -d $result_dir/06.enrich/loess ];then
  cp $result_dir/06.enrich/loess/loess.all.table ./03.loess/loess.all.table
  awk -v OFS='\t' '{$NF="";$(NF-1)="";$(NF-2)="";print $0}' ./03.loess/loess.all.table > ./03.loess/loess.all.table.xls
  rm -rf ./03.loess/loess.all.table
  cp $result_dir/06.enrich/loess/loess.region_stat.xls ./03.loess/loess.region_stat.xls
fi

# 04.enrich整理结果
for i in index loess;do
    mkdir -p ./04.enrich/$i
    cp -r $result_dir/06.enrich/$i/GO_result ./04.enrich/$i/
    cp -r $result_dir/06.enrich/$i/KEGG_result ./04.enrich/$i/
    cp $result_dir/06.enrich/$i.degfile ./04.enrich
    cp $result_dir/06.enrich/$i.genes_abstract.list ./04.enrich
done

## 生成BSA report
cd $workdir/
if [ ! -d ./report ];then
    mkdir -p ./report/Result
    mkdir -p ./report/file
fi

### 把data_release的结果放到report/Result
cp -r $workdir/data_release $workdir/report/Result
cp $workdir/data_info/group.list $workdir/report/file
cp $workdir/data_info/material_info.xls $workdir/report/file
cp $workdir/data_info/mix_info.xls $workdir/report/file
cp $workdir/data_info/project_info.xls $workdir/report/file
cp $workdir/data_info/stat_info.xls $workdir/report/file

if [ -e ./data_release/02.index/index.region_stat.xls ];then
  sed '1d' $workdir/data_release/02.index/index.region_stat.xls | sed 's/^/index-slid\t/g' > $workdir/report/file/index-slid.region_stat
fi

if [ -e ./data_release/03.loess/loess.region_stat.xls ];then
  sed '1d' $workdir/data_release/03.loess/loess.region_stat.xls | sed 's/^/index-loess\t/g' > $workdir/report/file/index-loess.region_stat
fi

cat $workdir/report/file/*.region_stat | sed '1i\Method\tRegion\tSNP_Number\tEffective_SNP\tINDEL_Number\tEffective_INDEL\tGene_Number\tNRID\tUniID\tKoID\tGOTERM\tEGGNOG\tPfamID'> $workdir/report/file/all.region_stat.xls

cp $workdir/data_release/01.vcf2table/snp_indel_gene.stat.xls $workdir/report/file/snp_indel_gene.stat.xls
cp $workdir/data_release/02.index/pop.index.index.png $workdir/report/file/slid_index.png
cp $workdir/data_release/03.loess/loess.index.png $workdir/report/file/loess_index.png

if [ -e $workdir/data_release/04.enrich/index/GO_result/index_GOenrichment.png ];then
  cp $workdir/data_release/04.enrich/index/GO_result/index_GOenrichment.png $workdir/report/file/index_go_enrich.png
fi

if [ -e $workdir/data_release/04.enrich/loess/GO_result/loess_GOenrichment.png ];then
  cp $workdir/data_release/04.enrich/loess/GO_result/loess_GOenrichment.png $workdir/report/file/loess_go_enrich.png
fi

if [ -e $workdir/data_release/04.enrich/index/KEGG_result/index_KEGGenrichment.png ];then
  cp $workdir/data_release/04.enrich/index/KEGG_result/index_KEGGenrichment.png $workdir/report/file/index_kegg_enrich.png
fi

if [ -e $workdir/data_release/04.enrich/loess/KEGG_result/loess_KEGGenrichment.png ];then
  cp $workdir/data_release/04.enrich/loess/KEGG_result/loess_KEGGenrichment.png $workdir/report/file/loess_kegg_enrich.png
fi

# 生辰报告
cp -r /mnt/lustre/users/sanger-dev/sg-users/yuan.xu/scripts/production_pipeline/rmarkdown_template_mutmap $workdir/report

## qc结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.raw.base.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.clean.base.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.raw.qual.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.clean.qual.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_qc_result.py --qc_stat $workdir/workflow_results/tmp/01.fastq_qc/qc.stat \
--group_file $workdir/report/file/group.list --raw_data $workdir/report/file/rawdata.xls --clean_data $workdir/report/file/cleandata.xls

## 比对结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workdir/workflow_results/tmp/03.mappingStat/$i.genome.coverage.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_mapping_result.py --all_summary_stats $workdir/workflow_results/tmp/03.mappingStat/all.summary.stats \
--group_file $workdir/report/file/group.list --bsa_align_stat $workdir/report/file/align_stat.xls --bsa_coverage_stat $workdir/report/file/coverage_sample.xls

## 变异检测结果整理
### SNP结果整理
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_snp_result.py --group_file $workdir/report/file/group.list \
--snp_stat $workdir/workflow_results/tmp/04.snpIndel/snp/snp.stat --snp_stat_all $workdir/workflow_results/tmp/04.snpIndel/snp/snp.stat.all \
--snp_stat_bsa $workdir/report/file/snp_stat.xls
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_indel_result.py --group_file $workdir/report/file/group.list \
--indel_stat $workdir/workflow_results/tmp/04.snpIndel/indel/indel.stat --indel_stat_bsa $workdir/report/file/indel_stat.xls
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_snp_anno_result.py --group_file $workdir/report/file/group.list \
--snp_anno_stat $workdir/workflow_results/tmp/04.snpIndel/snp/snp_anno.stat \
--snp_anno_stat_bsa $workdir/report/file/snp_anno.xls
cut -f 1,10,12,17-19,21 $workdir/report/file/snp_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/snp_anno_T.xls
cut -f 1,5-8 $workdir/report/file/snp_anno.xls > $workdir/report/file/snp_impact.xls
python3 $workdir/report/rmarkdown_template_mutmap/bin/process_indel_anno_result.py --group_file $workdir/report/file/group.list \
--indel_anno_stat $workdir/workflow_results/tmp/04.snpIndel/indel/indel_anno.stat \
--indel_anno_stat_bsa $workdir/report/file/indel_anno.xls
cut -f 1,16,17,19,27,29,30 $workdir/report/file/indel_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/indel_anno_T.xls
cut -f 1,6-9 $workdir/report/file/indel_anno.xls > $workdir/report/file/indel_impact.xls

cat $workdir/report/rmarkdown_template_mutmap/title.rmd $workdir/report/rmarkdown_template_mutmap/workflow.rmd $workdir/report/rmarkdown_template_mutmap/1.Fastq_qc.rmd \
$workdir/report/rmarkdown_template_mutmap/2.Align.rmd $workdir/report/rmarkdown_template_mutmap/3.SNP.rmd $workdir/report/rmarkdown_template_mutmap/4.Indel.rmd \
$workdir/report/rmarkdown_template_mutmap/7.亲本标记开发.rmd $workdir/report/rmarkdown_template_mutmap/index方法.rmd $workdir/report/rmarkdown_template_mutmap/候选区域的基因筛选.rmd \
$workdir/report/rmarkdown_template_mutmap/Go富集分析.rmd $workdir/report/rmarkdown_template_mutmap/KEGG富集分析.rmd $workdir/report/rmarkdown_template_mutmap/appendix.rmd > \
$workdir/report/rmarkdown_template_mutmap/bsa_all_mutmap.rmd

Rscript  $workdir/report/rmarkdown_template_mutmap/bin/rmarkdown --rmd  $workdir/report/rmarkdown_template_mutmap/bsa_all_mutmap.rmd \
--format html --outfile $workdir/report/bsa_all.html

rm -rf $workdir/report/rmarkdown_template
mv $workdir/report/bsa_all.html $workdir/report/Result

cp -r $workdir/report/rmarkdown_template_mutmap/readme_dir_mutmap $workdir/
python3 $workdir/readme_dir_mutmap/ReadMe_result.py --file $workdir/readme_dir_mutmap/BSA_all_readme.txt --typeinfo BSA_mutmap --outname BSA_mutmap
rm -rf $workdir/readme_dir_mutmap
mv $workdir/BSA_mutmap.ReadME.html $workdir/report/Result/data_release/

