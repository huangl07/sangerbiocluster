#!bin/bash
source  /mnt/lustre/users/sanger-dev/app/bioinfo/dna/new.rc
#是否运行过wgs  tdna结果路径 输入WGS结果路径/None

#task_id=$1
workdir=`pwd`
result_dir=$2
insert_name=$3
#workflow_results=`readNlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/`

#resultdir=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/tmp`
#filedir=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/published/data`

DNA_PATH=$(dirname $0)
echo "DNA_PATH $DNA_PATH"

cd $workdir

## 报告模板导入
if [ ! -d ./Report ];then
    mkdir -p ./Report/Result
    mkdir -p ./Report/info
fi

cp -r $DNA_PATH/rmd/ $workdir/Report
cp -r $DNA_PATH/readme/ $workdir/Report

#WGS报告
info_dir="$workdir/Report/info"
#cp $workdir/data_info/project.info $info_dir/project.info
echo $result_dir

if [ "$1" = "yes" ]; then
  workflow_results="$3/output/"
  echo $workflow_results
  resultdir="$3/output/tmp"
  filedir="$3/output/published/data"
  cp $resultdir/02.reference/project.info $info_dir/project.info
  cp $resultdir/02.reference/info.log $info_dir/info.log
  cp $resultdir/01.fastq_qc/qc.stat $info_dir
# cp $resultdir/04.snpIndel/snp/snp.stat.all $info_dir
# cp $resultdir/04.snpIndel/snp/snp.stat $info_dir
# cp $resultdir/04.snpIndel/snp/snp_anno.stat $info_dir
# cp $resultdir/04.snpIndel/indel/indel.stat $info_dir
  cp $resultdir/01.fastq_qc/json/* $info_dir
  cp $resultdir/03.mappingStat/all.summary.stats $info_dir
  cp $resultdir/03.mappingStat/*.png $info_dir
# awk '{print $1,$2,$3,$4,$5}' $resultdir/04.snpIndel/indel/indel.stat | grep -v "pop" |sed 's/ /\t/g'|sed 1d |sed '1i\Sample ID\tInsertion Number\tDeletion Number\tHeterozygosity Number\tHomozygosity Number'>$info_dir/indel_stat.xls
# cut -f1,16,17,20,27,29,30 $resultdir/04.snpIndel/indel/indel_anno.stat |sed 1d |sed '1i\Sample ID\tExon Loss Variant\tFrameshift Variant\tIntragenic Variant\tStart Lost\tStop Gained\tStop Lost'>$info_dir/indel_anno.xls
# cut -f1,6,7,8,9 $resultdir/04.snpIndel/indel/indel_anno.stat |sed 's/sampleID/Sample ID/g' >$info_dir/indel_anno_stat.xls
# cat $resultdir/09.anno/pop.stat.xls  |sed 's/#type/type/g' >$info_dir/gene_anno.xls
  cp -r $filedir/01.fastq_qc ./Report/Result
  cp -r $filedir/03.mappingStat ./Report/Result/
  mv ./Report/Result/03.mappingStat ./Report/Result/02.mappingStat
else
  echo "run no-wgs module"
#fq报告
  mkdir -p ./Report/Result/01.fastq_qc
  mkdir -p ./Report/Result/01.fastq_qc/json
  mkdir -p ./Report/Result/02.mappingStat
  cp $result_dir/fq/*.png ./Report/Result/01.fastq_qc/json
  cp $result_dir/fq/*.pdf ./Report/Result/01.fastq_qc/json
  cp $result_dir/fq/*.html ./Report/Result/01.fastq_qc/json
  cp $result_dir/fq/*.json ./Report/Result/01.fastq_qc/json
  cat $result_dir/fq/*.stat |grep -v '#' |cut -f1-13 |sort -r -u |sed '1i\#sampleID\trawreads\trawdata\trawq20\trawq30\trawGC\tada\tcleanreads\tcleandata\tcleanq20\tcleanq30\tcleanGC\tdup'> ./Report/Result/01.fastq_qc/qc.stat.txt
  #mapping_stat报告
  cp $result_dir/samtools_stat/*.png ./Report/Result/02.mappingStat/
  cp $result_dir/samtools_stat/*.pdf ./Report/Result/02.mappingStat/
  cat $result_dir/samtools_stat/*.summary.stat | head -1 > ./Report/Result/02.mappingStat/all.summary.stats.xls
  cat $result_dir/samtools_stat/*.summary.stat | grep -w -v "^Sample" | sort -u >> ./Report/Result/02.mappingStat/all.summary.stats.xls
  cp $result_dir/mosdepth_stat/*.pdf ./Report/Result/02.mappingStat/
  cp $result_dir/mosdepth_stat/*.png ./Report/Result/02.mappingStat/
  #整理info_dir
  cp ./Report/Result/02.mappingStat/all.summary.stats.xls $info_dir/all.summary.stats
  cp ./Report/Result/02.mappingStat/*.png $info_dir
  cp ./Report/Result/01.fastq_qc/qc.stat.txt $info_dir/qc.stat
  cp ./Report/Result/01.fastq_qc/json/* $info_dir
  cp $workdir/data/project.info $info_dir
  cp $workdir/data/info.log $info_dir
fi

#TDNA报告
mkdir -p ./Report/Result/03.tdna
mkdir -p ./Report/Result/03.tdna/tdnascan
mkdir -p ./Report/Result/03.tdna/aimhii
mkdir -p ./Report/Result/03.tdna/tdna_coverage
#mkdir -p ./Report/Result/03.tdna/tloc

#整理tdna_coverage结果文件目录
cp $result_dir/ref_cov/*.png  ./Report/Result/03.tdna/tdna_coverage/
for file in ./Report/Result/03.tdna/tdna_coverage/*.png; do
  mv "$file" "${file%.png}.tdna_coverage.png"
done

cp ./Report/Result/03.tdna/tdna_coverage/*.png $info_dir

#整理aimhii结果文件目录
if [ $# -lt 3 ];
then
    insert_name=$(awk 'FNR==NR{a[$0];next} {if(substr($0, 1, 1) == ">" && !($0 in a)) print $0}' <(grep ">" $result_dir/ref/ref.fa) <(grep ">" $result_dir/ref/insert_ref.fa) |cut -f2- -d">")
fi
#insert_name=`awk 'FNR==NR{a[$0];next} {if(substr($0, 1, 1) == ">" && !($0 in a)) print $0}' $result_dir/ref/ref.fa $result_dir/ref/insert_ref.fa|cut -f2- -d">"`
#insert_name=`grep -v -f $result_dir/ref/ref.fa $result_dir/ref/insert_ref.fa | grep ">" | cut -f2- -d">"`
echo "insert_name $insert_name"
python3 $DNA_PATH/tdna_result.py $result_dir/aimhii $info_dir/aimhii.stat.xls $insert_name

cp $info_dir/aimhii.stat.xls ./Report/Result/03.tdna/aimhii/
cp -r $result_dir/aimhii ./Report/Result/03.tdna/
rm -rf ./Report/Result/03.tdna/aimhii/*_temp
sed -i "/$insert_name/d" ./Report/Result/03.tdna/aimhii/*.csv
#整理tloc结果文件目录
#cp $result_dir/tloc/*/Result/Mosaic*.pdf ./Report/Result/03.tdna/tloc/
#cp $result_dir/tloc/*/TDNA/*mosaic_all.txt ./Report/Result/03.tdna/tloc/


#整理报告stat文件
rm -f $info_dir/tdna.stat.xls
touch $info_dir/tdna.stat.xls

echo "Sample	Chr	Position	SuppRead	TDNA_info	Orientation	Freq" >> $info_dir/tdna.stat.xls
for file in $result_dir/tdna/*/5.*.bed; do     filename=$(basename "$file");     sample_id=${filename/_insertion.bed/}; sample_id=${sample_id/5./}; awk -v id="$sample_id" 'NR>1{print id"\t"$0}' "$file" >> $info_dir/tdna.stat.xls; done
cat $info_dir/tdna.stat.xls
cp $info_dir/tdna.stat.xls ./Report/Result/03.tdna/tdnascan

cp -r $DNA_PATH/rmd/ $workdir/Report

cp -r $workdir/Report/info $workdir/Report/rmd/


bash $DNA_PATH/readme_generator.sh $workdir/Report/Result/
echo "try python"
cd $workdir/Report/Result/
python3 $workdir/Report/readme/readme.py -f WGS.result.info.txt -t WGS -n WGS
cd $workdir/Report/
cp $workdir/Report/readme/结果说明文档.txt $workdir/Report/Result/03.tdna
cp $workdir/Report/readme/01.readme.txt $workdir/Report/Result/01.fastq_qc/
mv $workdir/Report/Result/01.fastq_qc/01.readme.txt $workdir/Report/Result/01.fastq_qc/结果说明文档.txt
cp $workdir/Report/readme/02.readme.txt $workdir/Report/Result/02.mappingStat/
mv $workdir/Report/Result/02.mappingStat/02.readme.txt $workdir/Report/Result/02.mappingStat/结果说明文档.txt

mv $workdir/Report/Result/WGS.ReadME.html $workdir/Report/Result/结果目录索引及说明.html


Rscript $workdir/Report/rmd/rmarkdown.R --rmd $workdir/Report/rmd/report.rmd --format html --outfile $workdir/Report/tdna_all.html
Rscript $workdir/Report/rmd/rmarkdown.R --format pdf --outfile $workdir/Report/tdna_all.pdf --rmd $workdir/Report/rmd/report.rmd

mv $workdir/Report/tdna_all.html 插入位点Report.html
mv $workdir/Report/tdna_all.pdf 插入位点Report.pdf

rm -rf $workdir/Report/tdna_all.md
rm -rf $workdir/Report/tdna_all.tex
rm -rf $workdir/Report/rmd
rm -rf $workdir/Report/readme
rm -rf $workdir/Report/info
rm -rf $workdir/Report/Result/WGS.result.info.txt
