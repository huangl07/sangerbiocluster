#!bin/bash

workdir=`pwd`
data_realease_dir="$workdir/data_release"
result_dir="$workdir/result/"

cd $workdir
if [ ! -d $data_realease_dir ];then
    mkdir -p $data_realease_dir
fi

cd $data_realease_dir
mkdir 01.vcf2table 02.ridit 03.enrich

# 01.vcf2table结果整理
cp $result_dir/01.vcf2table/pop.table ./01.vcf2table/pop.table
cp $result_dir/01.vcf2table/snp_indel_gene.stat ./01.vcf2table/snp_indel_gene.stat.xls
cp $result_dir/01.vcf2table/pop.final.anno ./01.vcf2table/pop.final.anno
awk -v OFS='\t' '{$NF="";$(NF-1)="";$(NF-2)="";print $0}' ./01.vcf2table/pop.final.anno > ./01.vcf2table/pop.final.anno.xls
rm -rf ./01.vcf2table/pop.final.anno

# 02.ridit结果整理
cp $result_dir/02.ridit/pop.denoise.result ./02.ridit/pop.denoise.result.xls
cp $result_dir/02.ridit/loess.index.pdf ./02.ridit
cp $result_dir/02.ridit/loess.index.png ./02.ridit
cp $result_dir/02.ridit/loess.chr*.index.png ./02.ridit
cp $result_dir/02.ridit/loess.chr*.index.pdf ./02.ridit
if [ -d $result_dir/06.enrich/ridit ];then
  cp $result_dir/06.enrich/ridit/ridit.all.table ./02.ridit/ridit.all.table
  awk -v OFS='\t' '{$NF="";$(NF-1)="";$(NF-2)="";print $0}' ./02.ridit/ridit.all.table > ./02.ridit/ridit.all.table.xls
  rm -rf ./02.ridit/ridit.all.table
  cp $result_dir/06.enrich/ridit/ridit.region_stat.xls ./02.ridit/ridit.region_stat.xls
fi

# 03.enrich整理结果
for i in ridit;do
    mkdir -p ./03.enrich/$i
    cp -r $result_dir/06.enrich/$i/GO_result ./03.enrich/$i/
    cp -r $result_dir/06.enrich/$i/KEGG_result ./03.enrich/$i/
    cp $result_dir/06.enrich/$i.degfile ./03.enrich
    cp $result_dir/06.enrich/$i.genes_abstract.list ./03.enrich
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

if [ -e ./data_release/02.ridit/ridit.region_stat.xls ];then
  sed '1d' $workdir/data_release/02.ridit/ridit.region_stat.xls | sed 's/^/ridit\t/g' > $workdir/report/file/ridit.region_stat
fi

cat $workdir/report/file/*.region_stat | sed '1i\Method\tRegion\tSNP_Number\tEffective_SNP\tINDEL_Number\tEffective_INDEL\tGene_Number\tNRID\tUniID\tKoID\tGOTERM\tEGGNOG\tPfamID'> $workdir/report/file/all.region_stat.xls

cp $workdir/data_release/01.vcf2table/snp_indel_gene.stat.xls $workdir/report/file/snp_indel_gene.stat.xls
cp $workdir/data_release/02.ridit/loess.index.png $workdir/report/file/ridit.png

if [ -e $workdir/data_release/03.enrich/ridit/GO_result/ridit_GOenrichment.png ];then
  cp $workdir/data_release/03.enrich/ridit/GO_result/ridit_GOenrichment.png $workdir/report/file/ridit_go_enrich.png
fi

if [ -e $workdir/data_release/03.enrich/ridit/KEGG_result/ridit_KEGGenrichment.png ];then
  cp $workdir/data_release/03.enrich/ridit/KEGG_result/ridit_KEGGenrichment.png $workdir/report/file/ridit_kegg_enrich.png
fi

# 生辰报告
cp -r /mnt/lustre/users/sanger-dev/sg-users/yuan.xu/scripts/production_pipeline/rmarkdown_template_ridit $workdir/report

## qc结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.raw.base.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.clean.base.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.raw.qual.png $workdir/report/file
  cp $workdir/workflow_results/tmp/01.fastq_qc/json/$i.clean.qual.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_ridit/bin/process_qc_result.py --qc_stat $workdir/workflow_results/tmp/01.fastq_qc/qc.stat \
--group_file $workdir/report/file/group.list --raw_data $workdir/report/file/rawdata.xls --clean_data $workdir/report/file/cleandata.xls

## 比对结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workdir/workflow_results/tmp/03.mappingStat/$i.genome.coverage.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_ridit/bin/process_mapping_result.py --all_summary_stats $workdir/workflow_results/tmp/03.mappingStat/all.summary.stats \
--group_file $workdir/report/file/group.list --bsa_align_stat $workdir/report/file/align_stat.xls --bsa_coverage_stat $workdir/report/file/coverage_sample.xls

## 变异检测结果整理
### SNP结果整理
python3 $workdir/report/rmarkdown_template_ridit/bin/process_snp_result.py --group_file $workdir/report/file/group.list \
--snp_stat $workdir/workflow_results/tmp/04.snpIndel/snp/snp.stat --snp_stat_all $workdir/workflow_results/tmp/04.snpIndel/snp/snp.stat.all \
--snp_stat_bsa $workdir/report/file/snp_stat.xls
python3 $workdir/report/rmarkdown_template_ridit/bin/process_indel_result.py --group_file $workdir/report/file/group.list \
--indel_stat $workdir/workflow_results/tmp/04.snpIndel/indel/indel.stat --indel_stat_bsa $workdir/report/file/indel_stat.xls
python3 $workdir/report/rmarkdown_template_ridit/bin/process_snp_anno_result.py --group_file $workdir/report/file/group.list \
--snp_anno_stat $workdir/workflow_results/tmp/04.snpIndel/snp/snp_anno.stat \
--snp_anno_stat_bsa $workdir/report/file/snp_anno.xls
cut -f 1,10,12,17-19,21 $workdir/report/file/snp_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/snp_anno_T.xls
cut -f 1,5-8 $workdir/report/file/snp_anno.xls > $workdir/report/file/snp_impact.xls
python3 $workdir/report/rmarkdown_template_ridit/bin/process_indel_anno_result.py --group_file $workdir/report/file/group.list \
--indel_anno_stat $workdir/workflow_results/tmp/04.snpIndel/indel/indel_anno.stat \
--indel_anno_stat_bsa $workdir/report/file/indel_anno.xls
cut -f 1,16,17,19,27,29,30 $workdir/report/file/indel_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/indel_anno_T.xls
cut -f 1,6-9 $workdir/report/file/indel_anno.xls > $workdir/report/file/indel_impact.xls

cat $workdir/report/rmarkdown_template_ridit/title.rmd $workdir/report/rmarkdown_template_ridit/workflow.rmd $workdir/report/rmarkdown_template_ridit/1.Fastq_qc.rmd \
$workdir/report/rmarkdown_template_ridit/2.Align.rmd $workdir/report/rmarkdown_template_ridit/3.SNP.rmd $workdir/report/rmarkdown_template_ridit/4.Indel.rmd \
$workdir/report/rmarkdown_template_ridit/7.亲本标记开发.rmd $workdir/report/rmarkdown_template_ridit/ridit方法.rmd $workdir/report/rmarkdown_template_ridit/候选区域的基因筛选.rmd \
$workdir/report/rmarkdown_template_ridit/Go富集分析.rmd $workdir/report/rmarkdown_template_ridit/KEGG富集分析.rmd $workdir/report/rmarkdown_template_ridit/appendix.rmd > \
$workdir/report/rmarkdown_template_ridit/bsa_all_ridit.rmd

Rscript  $workdir/report/rmarkdown_template_ridit/bin/rmarkdown --rmd  $workdir/report/rmarkdown_template_ridit/bsa_all_ridit.rmd \
--format html --outfile $workdir/report/bsa_all.html

# 生成readme
cp -r $workdir/report/rmarkdown_template_mutmap/readme_dir_ridit $workdir/
echo "01.vcf2table	folder_1	变异检测分析结果
01.vcf2table_readme.txt	file	readme文件
pop.table	file	变异检测统计列表
pop.final.anno.xls	file	变异检测带基因注释表统计列表
snp_indel_gene.stat.xls	file	变异检测统计列表统计
02.ridit	folder_1	ridit方法获取的定位结果文件夹
02.ridit_readme.txt	file	readme文件
pop.denoise.result.xls	file	ridit方法结果文件
ridit.region_stat.xls	file	ridit方法定位区域统计表
ridit.all.table.xls	file	ridit方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "loess.$i.index.png	file	不同染色体上ridit方法的曼哈顿图（png格式）
loess.$i.index.pdf	file	不同染色体上ridit方法的曼哈顿图（pdf格式）" >> BSA.info.result
done

echo "loess.index.pdf	file	ridit方法的整体曼哈顿图（pdf格式）
loess.index.png	file	ridit方法的整体曼哈顿图（png格式）
03.enrich	folder_1	候选基因GO和KEGG富集分析
03.enrich_readme.txt	file	readme文件
ridit	folder_2	候选基因GO和KEGG富集分析（index方法）
ridit.genes_abstract.list	file	定位区域基因bed文件
ridit.degfile	file	定位区域基因列表
GO_result	folder_3	候选基因GO富集结果
ridit_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
ridit_GOenrichment.xls	file	候选基因GO富集结果表
ridit_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
ridit_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
ridit_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
ridit_KEGGenrichment.xls	file	候选基因KEGG富集结果表
ridit_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
ridit_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）" >> BSA.info.result

cp -r $workdir/report/rmarkdown_template_ridit/readme_dir_ridit $workdir/
python3 $workdir/readme_dir_ridit/ReadMe_result.py --file BSA.info.result --typeinfo BSA_ridit --outname BSA_ridit
mv $workdir/BSA_ridit.ReadME.html $workdir/report/Result/data_release/

rm -rf $workdir/report/rmarkdown_template_ridit
rm -rf $workdir/readme_dir_ridit
rm -rf BSA.info.result
mv $workdir/report/bsa_all.html $workdir/report/Result
