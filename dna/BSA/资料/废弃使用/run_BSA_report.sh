#!bin/bash
###
 # @Author: error: git config user.name && git config user.email & please set dead value or install git
 # @Date: 2022-12-06 13:45:42
 # @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @LastEditTime: 2023-02-08 14:13:11
 # @FilePath: /ZQH_5612/mnt/ilustre/users/jing.wang/work/FX2022112400238/test/run_report.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

task_id=$1
workdir=`pwd`
data_realease_dir="$workdir/data_release"
result_dir="$workdir/result/"
workflow_results=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/`
echo $workflow_results


cd $workdir
if [ ! -d $data_realease_dir ];then
    mkdir -p $data_realease_dir
fi

## 报告模板导入
if [ ! -d ./report ];then
    mkdir -p ./report/Result
fi
ln -s $workdir/data/project.txt $workdir/report
ln -s $workdir/data/project.info.xls $workdir/report/Result
ln -s $workdir/data/result.info.xls $workdir/report/Result
ln -s $workdir/data/chr.list $workdir/report/Result
ln -s $workdir/workflow_results $workdir/report/Result
cp -r $workdir/result/ $workdir/report/Result/Result_bsa
cat $data_realease_dir/03.index-slid/gene.total.xls |grep "@" |cut -f 1-5 |sed 's/@//g' |awk '{print$1"_"$2"_"$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' >$workdir/report/Result/Result_bsa/02.index-slid/index.gene.total
if [ -s $data_realease_dir/02.ED-slid/gene.total.xls ];then
cat $data_realease_dir/02.ED-slid/gene.total.xls |grep "@" |cut -f 1-5 |sed 's/@//g' |awk '{print$1"_"$2"_"$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' >$workdir/report/Result/Result_bsa/02.index-slid/ED.gene.total
fi
cp -r /mnt/lustre/users/sanger-dev/sg-users/yuan.xu/scripts/production_pipeline/rmarkdown_template_pdf $workdir/report
cp $workflow_results/tmp/02.reference/project.info $workdir/report/file
cp $workflow_results/tmp/02.reference/info.log $workdir/report/file

cd $data_realease_dir
mkdir 01.vcf2table  02.index  03.loess  04.ED  05.Gprime  06.DeepBSA_DL 07.DeepBSA_K 08.enrich

# 01.vcf2table结果整理
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $result_dir/01.vcf2table/pop.index | cut -f 3-5 > tmp2.txt
cat $result_dir/01.vcf2table/pop.table |cut -f 1-5 > tmp1.txt
awk -F "\t" '{for (i=6;i<=NF;i++)printf("%s\t", $i);print ""}' $result_dir/01.vcf2table/pop.table > tmp3.txt
paste tmp1.txt tmp2.txt tmp3.txt > ./01.vcf2table/pop.table
rm tmp1.txt tmp2.txt tmp3.txt
cp $result_dir/01.vcf2table/snp_indel_gene.stat ./01.vcf2table/snp_indel_gene.stat.xls
sed -i '1s/.*/CHROM\tSNP Number\tEffective SNP\tInDel Number\tEffective InDel\tGene Number\tEffective Gene/' ./01.vcf2table/snp_indel_gene.stat.xls
cp $result_dir/01.vcf2table/pop.final.anno ./01.vcf2table/pop.final.anno
Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile ./01.vcf2table/pop.final.anno --outfile ./01.vcf2table/pop.final.anno.xls
rm -rf ./01.vcf2table/pop.final.anno

# 02.index-slid结果整理
awk '{print $1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $result_dir/02.index-slid/pop.bootstrap.result > ./02.index/pop.index.bootstrap.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/02.index-slid/pop.index.$i.index.png ./02.index/index.$i.png
  cp $result_dir/02.index-slid/pop.index.$i.index.pdf ./02.index/index.$i.pdf
done
cp $result_dir/02.index-slid/pop.index.index.png ./02.index/index.png
cp $result_dir/02.index-slid/pop.index.index.pdf ./02.index/index.pdf
if [ -d $result_dir/06.enrich/index ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/index/index.all.table --outfile ./02.index/index.all.table.xls
  cp $result_dir/06.enrich/index/index.region_stat.xls ./02.index/index.region_stat.xls
fi

# 03.loess结果整理
awk '{print $1"\t"$2"\t"$8"\t"$9"\t"$10}' $result_dir/05.loess/pop.bootstrap.result >./03.loess/pop.loess.bootstrap.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/05.loess/loess.$i.index.png ./03.loess/loess.$i.png
  cp $result_dir/05.loess/loess.$i.index.pdf ./03.loess/loess.$i.pdf
done
cp $result_dir/05.loess/loess.index.png ./03.loess/loess.png
cp $result_dir/05.loess/loess.index.pdf ./03.loess/loess.pdf
if [ -d $result_dir/06.enrich/loess ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/loess/loess.all.table --outfile ./03.loess/loess.all.table.xls
  cp $result_dir/06.enrich/loess/loess.region_stat.xls ./03.loess/loess.region_stat.xls
fi

# 04.ED-slid结果整理
awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$8}' $result_dir/04.ED-slid/pop.sliding.result > ./04.ED/pop.ED.sliding.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/04.ED-slid/pop.ED.$i.index.png ./04.ED/ED.$i.png
  cp $result_dir/04.ED-slid/pop.ED.$i.index.pdf ./04.ED/ED.$i.pdf
done
cp $result_dir/04.ED-slid/pop.ED.index.png ./04.ED/ED.png
cp $result_dir/04.ED-slid/pop.ED.index.pdf ./04.ED/ED.pdf
if [ -d $result_dir/06.enrich/ED ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/ED/ED.all.table --outfile ./04.ED/ED.all.table.xls
  cp $result_dir/06.enrich/ED/ED.region_stat.xls ./04.ED/ED.region_stat.xls
fi

# 05.Gprime结果整理
awk '{print $1"\t"$2"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17}' $result_dir/03.Gprime/pop.Gprime > ./05.Gprime/pop.Gprime.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/03.Gprime/Gprime.$i.index.png ./05.Gprime/Gprime.$i.png
  cp $result_dir/03.Gprime/Gprime.$i.index.pdf ./05.Gprime/Gprime.$i.pdf
done
cp $result_dir/03.Gprime/Gprime.index.png ./05.Gprime/Gprime.png
cp $result_dir/03.Gprime/Gprime.index.pdf ./05.Gprime/Gprime.pdf
if [ -d $result_dir/06.enrich/Gprime ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/Gprime/Gprime.all.table --outfile ./05.Gprime/Gprime.all.table.xls
  cp $result_dir/06.enrich/Gprime/Gprime.region_stat.xls ./05.Gprime/Gprime.region_stat.xls
fi

# 06.DeepBSA_DL结果整理
cp $result_dir/08.DeepBSA/DL_values.txt ./06.DeepBSA_DL/pop.DeepBSA_DL.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/08.DeepBSA/pop.deepBSA_DL.$i.index.png ./06.DeepBSA_DL/DeepBSA_DL.$i.png
  cp $result_dir/08.DeepBSA/pop.deepBSA_DL.$i.index.pdf ./06.DeepBSA_DL/DeepBSA_DL.$i.pdf
done
cp $result_dir/08.DeepBSA/pop.deepBSA_DL.index.png ./06.DeepBSA_DL/DeepBSA_DL.png
cp $result_dir/08.DeepBSA/pop.deepBSA_DL.index.pdf ./06.DeepBSA_DL/DeepBSA_DL.pdf
if [ -d $result_dir/06.enrich/DeepBSA_DL ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/DeepBSA_DL/DeepBSA_DL.all.table --outfile ./06.DeepBSA_DL/DeepBSA_DL.all.table.xls
  cp $result_dir/06.enrich/DeepBSA_DL/DeepBSA_DL.region_stat.xls ./06.DeepBSA_DL/DeepBSA_DL.region_stat.xls
fi

# 07.DeepBSA_K结果整理
cp $result_dir/08.DeepBSA/K_values.txt ./07.DeepBSA_K/pop.DeepBSA_K.result.xls
for i in $(cut -f1 $workdir/data_info/chr.list);do
  cp $result_dir/08.DeepBSA/pop.deepBSA_K.$i.index.png ./07.DeepBSA_K/DeepBSA_K.$i.png
  cp $result_dir/08.DeepBSA/pop.deepBSA_K.$i.index.pdf ./07.DeepBSA_K/DeepBSA_K.$i.pdf
done
cp $result_dir/08.DeepBSA/pop.deepBSA_K.index.png ./07.DeepBSA_K/DeepBSA_K.png
cp $result_dir/08.DeepBSA/pop.deepBSA_K.index.pdf ./07.DeepBSA_K/DeepBSA_K.pdf
if [ -d $result_dir/06.enrich/DeepBSA_K ];then
  Rscript $workdir/report/rmarkdown_template_pdf/bin/process_final_anno.r --infile $result_dir/06.enrich/DeepBSA_K/DeepBSA_K.all.table --outfile ./07.DeepBSA_K/DeepBSA_K.all.table.xls
  cp $result_dir/06.enrich/DeepBSA_K/DeepBSA_K.region_stat.xls ./07.DeepBSA_K/DeepBSA_K.region_stat.xls
fi

# 08.enrich整理结果
for i in ED Gprime index loess DeepBSA_DL DeepBSA_K;do
    mkdir -p ./08.enrich/$i
    cp -r $result_dir/06.enrich/$i/GO_result ./08.enrich/$i/
    cp -r $result_dir/06.enrich/$i/KEGG_result ./08.enrich/$i/
    cp $result_dir/06.enrich/$i.degfile ./08.enrich/$i/$i.degfile.xls
    cp $result_dir/06.enrich/$i.transcript.degfile ./08.enrich/$i/$i.transcript.degfile.xls
    cp $result_dir/06.enrich/$i.genes_abstract.list ./08.enrich/$i/$i.genes_abstract.list.xls
done

### 把data_release的结果放到report/Result
cd $workdir
cp -r $workdir/data_release $workdir/report/Result
cp $workdir/data_info/group.list $workdir/report/file
cp $workdir/data_info/chr.list $workdir/report/file

if [ -e ./data_release/02.index/index.region_stat.xls ];then
  sed '1d' $workdir/data_release/02.index/index.region_stat.xls | sed 's/^/index-slid\t/g' > $workdir/report/file/index.region_stat
fi

if [ -e ./data_release/03.loess/loess.region_stat.xls ];then
  sed '1d' $workdir/data_release/03.loess/loess.region_stat.xls | sed 's/^/index-loess\t/g' > $workdir/report/file/loess.region_stat
fi

if [ -e ./data_release/04.ED/ED.region_stat.xls ];then
  sed '1d' $workdir/data_release/04.ED/ED.region_stat.xls | sed 's/^/Euclidean\t/g' > $workdir/report/file/ED.region_stat
fi

if [ -e ./data_release/05.Gprime/Gprime.region_stat.xls ];then
  sed '1d' $workdir/data_release/05.Gprime/Gprime.region_stat.xls | sed 's/^/Gprime\t/g' > $workdir/report/file/Gprime.region_stat
fi

if [ -e ./data_release/06.DeepBSA_DL/DeepBSA_DL.region_stat.xls ];then
  sed '1d' $workdir/data_release/06.DeepBSA_DL/DeepBSA_DL.region_stat.xls | sed 's/^/DeepBSA_DL\t/g' > $workdir/report/file/DeepBSA_DL.region_stat
fi

if [ -e ./data_release/07.DeepBSA_K/DeepBSA_K.region_stat.xls ];then
  sed '1d' $workdir/data_release/07.DeepBSA_K/DeepBSA_K.region_stat.xls | sed 's/^/DeepBSA_K\t/g' > $workdir/report/file/DeepBSA_K.region_stat
fi

cat $workdir/report/file/*.region_stat | sed '1i\Method\tRegion\tSNP_Number\tEffective_SNP\tInDel_Number\tEffective_InDel\tGene_Number\tNRID\tUniID\tKoID\tGOTERM\tEGGNOG\tPfamID'> $workdir/report/file/all.region_stat.xls

cp $workdir/data_release/01.vcf2table/snp_indel_gene.stat.xls $workdir/report/file/snp_indel_gene.stat.xls
cp $workdir/data_release/02.index/index.png $workdir/report/file/
cp $workdir/data_release/03.loess/loess.png $workdir/report/file/
cp $workdir/data_release/04.ED/ED.png $workdir/report/file/
cp $workdir/data_release/05.Gprime/Gprime.png $workdir/report/file/
cp $workdir/data_release/06.DeepBSA_DL/DeepBSA_DL.png $workdir/report/file/
cp $workdir/data_release/07.DeepBSA_K/DeepBSA_K.png $workdir/report/file/

if [ -e $workdir/data_release/08.enrich/index/GO_result/index_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/index/GO_result/index_GOenrichment.png $workdir/report/file/index_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/loess/GO_result/loess_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/loess/GO_result/loess_GOenrichment.png $workdir/report/file/loess_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/ED/GO_result/ED_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/ED/GO_result/ED_GOenrichment.png $workdir/report/file/ED_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/Gprime/GO_result/Gprime_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/Gprime/GO_result/Gprime_GOenrichment.png $workdir/report/file/Gprime_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/DeepBSA_DL/GO_result/DeepBSA_DL_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/DeepBSA_DL/GO_result/DeepBSA_DL_GOenrichment.png $workdir/report/file/DeepBSA_DL_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/DeepBSA_K/GO_result/DeepBSA_K_GOenrichment.png ];then
  cp $workdir/data_release/08.enrich/DeepBSA_K/GO_result/DeepBSA_K_GOenrichment.png $workdir/report/file/DeepBSA_K_go_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/index/KEGG_result/index_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/index/KEGG_result/index_KEGGenrichment.png $workdir/report/file/index_kegg_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/loess/KEGG_result/loess_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/loess/KEGG_result/loess_KEGGenrichment.png $workdir/report/file/loess_kegg_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/ED/KEGG_result/ED_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/ED/KEGG_result/ED_KEGGenrichment.png $workdir/report/file/ED_kegg_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/Gprime/KEGG_result/Gprime_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/Gprime/KEGG_result/Gprime_KEGGenrichment.png $workdir/report/file/Gprime_kegg_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/DeepBSA_DL/KEGG_result/DeepBSA_DL_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/DeepBSA_DL/KEGG_result/DeepBSA_DL_KEGGenrichment.png $workdir/report/file/DeepBSA_DL_kegg_enrich.png
fi

if [ -e $workdir/data_release/08.enrich/DeepBSA_K/KEGG_result/DeepBSA_K_KEGGenrichment.png ];then
  cp $workdir/data_release/08.enrich/DeepBSA_K/KEGG_result/DeepBSA_K_KEGGenrichment.png $workdir/report/file/DeepBSA_K_kegg_enrich.png
fi

cp $workdir/data_release/08.enrich/*/*.degfile.xls $workdir/report/file

## qc结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workflow_results/tmp/01.fastq_qc/json/$i.raw.base.png $workdir/report/file
  cp $workflow_results/tmp/01.fastq_qc/json/$i.clean.base.png $workdir/report/file
  cp $workflow_results/tmp/01.fastq_qc/json/$i.raw.qual.png $workdir/report/file
  cp $workflow_results/tmp/01.fastq_qc/json/$i.clean.qual.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_pdf/bin/process_qc_result.py --qc_stat $workflow_results/tmp/01.fastq_qc/qc.stat \
--group_file $workdir/report/file/group.list --raw_data $workdir/report/file/rawdata.xls --clean_data $workdir/report/file/cleandata.xls

## 比对结果整理
for i in `cut -f 1 $workdir/report/file/group.list`
do
  cp $workflow_results/tmp/03.mappingStat/$i.genome.coverage.png $workdir/report/file
  cp $workflow_results/tmp/03.mappingStat/$i.insert.png $workdir/report/file
  cp $workflow_results/tmp/03.mappingStat/$i.depth.png $workdir/report/file
done
python3 $workdir/report/rmarkdown_template_pdf/bin/process_mapping_result.py --all_summary_stats $workflow_results/tmp/03.mappingStat/all.summary.stats \
--group_file $workdir/report/file/group.list --bsa_align_stat $workdir/report/file/align_stat.xls

## 变异检测结果整理
### SNP结果整理
python3 $workdir/report/rmarkdown_template_pdf/bin/process_snp_result.py --group_file $workdir/report/file/group.list \
--snp_stat $workflow_results/tmp/04.snpIndel/snp/snp.stat --snp_stat_bsa $workdir/report/file/snp_stat.xls
python3 $workdir/report/rmarkdown_template_pdf/bin/process_indel_result.py --group_file $workdir/report/file/group.list \
--indel_stat $workflow_results/tmp/04.snpIndel/indel/indel.stat --indel_stat_bsa $workdir/report/file/indel_stat.xls
python3 $workdir/report/rmarkdown_template_pdf/bin/process_snp_anno_result.py --group_file $workdir/report/file/group.list \
--snp_anno_stat $workflow_results/tmp/04.snpIndel/snp/snp_anno.stat \
--snp_anno_stat_bsa $workdir/report/file/snp_anno.xls
cut -f 1,10,12,17-19,21 $workdir/report/file/snp_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/snp_anno_T.xls
cut -f 1,5-8 $workdir/report/file/snp_anno.xls > $workdir/report/file/snp_impact.xls
python3 $workdir/report/rmarkdown_template_pdf/bin/process_indel_anno_result.py --group_file $workdir/report/file/group.list \
--indel_anno_stat $workflow_results/tmp/04.snpIndel/indel/indel_anno.stat \
--indel_anno_stat_bsa $workdir/report/file/indel_anno.xls
cut -f 1,16,17,19,27,29,30 $workdir/report/file/indel_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
| sed 's/[ \t]*$//g' | sed 's/ /\t/g' > $workdir/report/file/indel_anno_T.xls
cut -f 1,6-9 $workdir/report/file/indel_anno.xls > $workdir/report/file/indel_impact.xls
cp $workdir/report/file/snp_anno.xls $workdir/report/Result/data_release/01.vcf2table
cp $workdir/report/file/indel_anno.xls $workdir/report/Result/data_release/01.vcf2table

Rscript $workdir/report/rmarkdown_template_pdf/rmarkdown.r --rmd  $workdir/report/rmarkdown_template_pdf/report.rmd \
--format html --outfile $workdir/report/bsa_all.html
Rscript $workdir/report/rmarkdown_template_pdf/rmarkdown.r --rmd  $workdir/report/rmarkdown_template_pdf/report.rmd \
--format pdf --outfile $workdir/report/bsa_all.pdf

echo "01.vcf2table	folder_1	变异检测分析结果
01.vcf2table_readme.txt	file	readme文件
pop.table	file	变异检测统计列表
pop.final.anno.xls	file	变异检测带基因注释表统计列表
snp_indel_gene.stat.xls	file	变异检测统计表
snp_anno.xls	file	SNP功效信息统计表
indel_anno.xls	file	InDel功效信息统计表
02.index	folder_1	index方法获取的定位结果文件夹
02.index_readme.txt	file	readme文件
pop.index.bootstrap.result.xls	file	index方法结果文件
index.region_stat.xls	file	index方法定位区域统计表
index.all.table.xls	file	index方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "index.$i.png	file	不同染色体上index方法的曼哈顿图（png格式）
index.$i.pdf	file	不同染色体上index方法的曼哈顿图（pdf格式）" >> BSA.info.result
done

echo "index.png	file	index方法的整体曼哈顿图（png格式）
index.pdf	file	index方法的整体曼哈顿图（pdf格式）
03.loess	folder_1	loess方法获取的定位结果文件夹
03.loess_readme.txt	file	readme文件
pop.loess.bootstrap.result.xls	file	loess方法结果文件
loess.region_stat.xls	file	loess方法定位区域统计表
loess.all.table.xls	file	loess方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "loess.$i.pdf	file	不同染色体上loess方法的曼哈顿图（pdf格式）
loess.$i.png	file	不同染色体上loess方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "loess.pdf	file	loess方法的整体曼哈顿图（pdf格式）
loess.png	file	loess方法的整体曼哈顿图（png格式）
04.ED	folder_1	ED方法获取的定位结果文件夹
04.ED_readme.txt	file	readme文件
pop.ED.sliding.result.xls	file	ED方法统计结果文件
ED.region_stat.xls	file	ED方法定位区域位点和基因详情表
ED.all.table.xls	file	ED方法的整体曼哈顿图（png格式）" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "ED.$i.pdf	file	不同染色体上ED方法的曼哈顿图（pdf格式）
ED.$i.png	file	不同染色体上ED方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "ED.pdf	file	ED方法的整体曼哈顿图（pdf格式）
ED.png	file	ED方法的整体曼哈顿图（png格式）
05.Gprime	folder_1	Gprime方法获取的定位结果文件夹
05.Gprime_readme.txt	file	readme文件
pop.Gprime.result.xls	file	Gprime方法统计结果文件
Gprime.region_stat.xls	file	Gprime方法定位区域统计表
Gprime.all.table.xls	file	Gprime方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "Gprime.$i.pdf	file	不同染色体上Gprime方法的曼哈顿图（pdf格式）
Gprime.$i.pdf	file	不同染色体上Gprime方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "Gprime.pdf	file	Gprime方法的整体曼哈顿图（pdf格式）
Gprime.png	file	Gprime方法的整体曼哈顿图（png格式）
06.DeepBSA_DL	folder_1	DeepBSA_DL方法获取的定位结果文件夹
06.DeepBSA_DL_readme.txt	file	readme文件
pop.DeepBSA_DL.result.xls	file	DeepBSA_DL方法统计结果文件
DeepBSA_DL.region_stat.xls	file	DeepBSA_DL方法定位区域统计表
DeepBSA_DL.all.table.xls	file	DeepBSA_DL方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "DeepBSA_DL.$i.pdf	file	不同染色体上DeepBSA_DL方法的曼哈顿图（pdf格式）
DeepBSA_DL.$i.pdf	file	不同染色体上DeepBSA_DL方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "DeepBSA_DL.pdf	file	DeepBSA_DL方法的整体曼哈顿图（pdf格式）
DeepBSA_DL.png	file	DeepBSA_DL方法的整体曼哈顿图（png格式）
07.DeepBSA_K	folder_1	DeepBSA_DL方法获取的定位结果文件夹
07.DeepBSA_K_readme.txt	file	readme文件
pop.DeepBSA_K.result.xls	file	DeepBSA_K方法统计结果文件
DeepBSA_K.region_stat.xls	file	DeepBSA_K方法定位区域统计表
DeepBSA_K.all.table.xls	file	DeepBSA_K方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "DeepBSA_K.$i.pdf	file	不同染色体上DeepBSA_K方法的曼哈顿图（pdf格式）
DeepBSA_K.$i.pdf	file	不同染色体上DeepBSA_K方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "DeepBSA_K.pdf	file	DeepBSA_K方法的整体曼哈顿图（pdf格式）
DeepBSA_K.png	file	DeepBSA_K方法的整体曼哈顿图（png格式）
08.enrich	folder_1	候选基因GO和KEGG富集分析
08.enrich_readme.txt	file	readme文件
index	folder_2	候选基因GO和KEGG富集分析（index方法）
index.degfile.xls	file	差异基因列表
index.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
index_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
index_GOenrichment.xls	file	候选基因GO富集结果表
index_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
index_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
index_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
index_KEGGenrichment.xls	file	候选基因KEGG富集结果表
index_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
index_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）
loess	folder_2	候选基因GO和KEGG富集分析（loess方法）
loess.degfile.xls	file	差异基因列表
loess.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
loess_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
loess_GOenrichment.xls	file	候选基因GO富集结果表
loess_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
loess_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
loess_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
loess_KEGGenrichment.xls	file	候选基因KEGG富集结果表
loess_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
loess_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）
ED	folder_2	候选基因GO和KEGG富集分析（ED方法）
ED.degfile.xls	file	差异基因列表
ED.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
ED_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
ED_GOenrichment.xls	file	候选基因GO富集结果表
ED_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
ED_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
ED_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
ED_KEGGenrichment.xls	file	候选基因KEGG富集结果表
ED_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
ED_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）
Gprime	folder_2	候选基因GO和KEGG富集分析（Gprime方法）
Gprime.degfile.xls	file	差异基因列表
Gprime.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
Gprime_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
Gprime_GOenrichment.xls	file	候选基因GO富集结果表
Gprime_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
Gprime_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
Gprime_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
Gprime_KEGGenrichment.xls	file	候选基因KEGG富集结果表
Gprime_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
Gprime_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）
DeepBSA_DL	folder_2	候选基因GO和KEGG富集分析（loess方法）
DeepBSA_DL.degfile.xls	file	差异基因列表
DeepBSA_DL.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
DeepBSA_DL_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
DeepBSA_DL_GOenrichment.xls	file	候选基因GO富集结果表
DeepBSA_DL_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
DeepBSA_DL_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
DeepBSA_DL_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
DeepBSA_DL_KEGGenrichment.xls	file	候选基因KEGG富集结果表
DeepBSA_DL_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
DeepBSA_DL_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）
DeepBSA_K	folder_2	候选基因GO和KEGG富集分析（loess方法）
DeepBSA_K.degfile.xls	file	差异基因列表
DeepBSA_K.genes_abstract.list.xls	file	定位区域基因bed文件
GO_result	folder_3	候选基因GO富集结果
DeepBSA_K_GOenrichment_0.05.xls	file	候选基因GO显著富集结果表
DeepBSA_K_GOenrichment.xls	file	候选基因GO富集结果表
DeepBSA_K_GOenrichment.pdf	file	候选区域内基因GO分类统计图（pdf格式）
DeepBSA_K_GOenrichment.png	file	候选区域内基因GO分类统计图（png格式）
KEGG_result	folder_3	候选基因KEGG富集结果
DeepBSA_K_KEGGenrichment_0.05.xls	file	候选基因KEGG显著富集结果表
DeepBSA_K_KEGGenrichment.xls	file	候选基因KEGG富集结果表
DeepBSA_K_KEGGenrichment.pdf	file	候选区域内基因的KEGG富集统计图（pdf格式）
DeepBSA_K_KEGGenrichment.png	file	候选区域内基因的KEGG富集统计图（png格式）" >> BSA.info.result

cp -r $workdir/report/rmarkdown_template_pdf/readme_dir_qtl $workdir/
python3 $workdir/readme_dir_qtl/ReadMe_result.py --file BSA.info.result --typeinfo BSA_new --outname BSA

rm -rf BSA.info.result
rm -rf $workdir/report/rmarkdown_template_pdf
rm -rf $workdir/readme_dir_qtl
rm -rf $workdir/data_release
mv $workdir/report/bsa_all.html $workdir/report/Result
mv $workdir/report/bsa_all.pdf $workdir/report/Result
mv $workdir/BSA.ReadME.html $workdir/report/Result/data_release
