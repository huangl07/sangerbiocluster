#$1 work_dir

output_data=$1
touch $output_data/noref.result.info.txt

cd $output_data
##############01
echo -e "01.stat\\tfolder_1\\t统计结果目录" >> $output_data/noref.result.info.txt
echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/noref.result.info.txt
echo -e "qc.stat.xls\\tfile\\t原始数据及质控统计表" >> $output_data/noref.result.info.txt
echo -e "sample_tags.stat.xls\\tfile\\t样品tags统计表" >> $output_data/noref.result.info.txt
echo -e "fig\\tfolder_2\\t样品质控结果汇总" >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.clean.base.pdf file $j样本指控后碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.clean.base.png file $j样本指控后碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.clean.qual.pdf file $j样本指控后碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.clean.qual.png file $j样本指控后碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.raw.base.pdf file $j样本原始碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.raw.base.png file $j样本原始碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.raw.qual.pdf file $j样本原始碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.raw.qual.png file $j样本原始碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt
for j in `ls 01.stat/fig |cut -d "." -f1 |sort -u`;do echo $j.json file $j样本质控结果; done |sed 's/ /\t/g' >> $output_data/noref.result.info.txt

#####02
echo -e "02.variants\\tfolder_1\\t变异检测结果汇总" >> $output_data/noref.result.info.txt
echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/noref.result.info.txt
echo -e "ref.loci.fa\\tfile\\t参考tags序列文件" >> $output_data/noref.result.info.txt
echo -e "ref.loci.fa.fai\\tfile\\t参考tags序列的索引文件" >> $output_data/noref.result.info.txt
echo -e "populations.final.vcf.gz\\tfile\\t样品变异检测结果文件" >> $output_data/noref.result.info.txt
echo -e "populations.final.vcf.gz.tbi\\tfile\\t样品变异检测结果的索引文件" >> $output_data/noref.result.info.txt
echo -e "snp.stat.xls\\tfile\\t变异基因型文件单样品snp统计表" >> $output_data/noref.result.info.txt
echo -e "snp.stat.all.xls\\tfile\\t群体变异基因型snp统计表" >> $output_data/noref.result.info.txt
