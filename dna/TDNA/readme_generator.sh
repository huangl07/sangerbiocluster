#$1 work_dir

output_data=$1

touch $output_data/WGS.result.info.txt
cd $output_data/
##############01
echo -e "01.fastq_qc\\tfolder_1\\t数据质控结果目录" >> $output_data/WGS.result.info.txt
echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
echo -e "qc.stat\\tfile\\t原始数据及质控统计表" >> $output_data/WGS.result.info.txt
echo -e "json\\tfolder_2\\t样品质控结果汇总" >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.clean.base.pdf file $j样本指控后碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.clean.base.png file $j样本指控后碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.clean.qual.pdf file $j样本指控后碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.clean.qual.png file $j样本指控后碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.raw.base.pdf file $j样本原始碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.raw.base.png file $j样本原始碱基组成分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.raw.qual.pdf file $j样本原始碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.raw.qual.png file $j样本原始碱基错误率分布例图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.json file $j样本质控结果; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt

#####02
##echo -e "02.reference\\tfolder_1\\t参考基因组文件结果汇总" >> $output_data/WGS.result.info.txt
##echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
##echo -e "ref.fa\\tfile\\t参考基因组文件" >> $output_data/WGS.result.info.txt
##echo -e "ref.gtf\\tfile\\t参考基因组注释文件" >> $output_data/WGS.result.info.txt
##echo -e "ref.cds.fa\\tfile\\t参考基因组（编码区）文件" >> $output_data/WGS.result.info.txt
##echo -e "ref.protein.fa\\tfile\\t参考基因组（蛋白）文件" >> $output_data/WGS.result.info.txt
##echo -e "ref.genome.summary.xls\\tfile\\t参考基因组序列统计表" >> $output_data/WGS.result.info.txt

########03
echo -e "02.mappingStat\\tfolder_1\\t比对统计结果汇总" >> $output_data/WGS.result.info.txt
echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
echo -e "all.summary.stats\\tfile\\t所有样品比对统计结果汇总表" >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.depth.pdf file $j样本基因组比对测序深度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.depth.png file $j样本基因组比对测序深度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.coverage.pdf file $j样本基因组比对覆盖度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.coverage.png file $j样本基因组比对覆盖度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.insert.pdf file $j样本基因组比对插入片段长度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
for j in `ls 02.mappingStat |cut -d "." -f1 |sort -u`;do echo $j.insert.png file $j样本基因组比对插入片段长度结果目录; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt

############04
#echo -e "04.snp_indel\\tfolder_1\\tsnp/indel变异检测统计结果汇总" >> $output_data/WGS.result.info.txt
#echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#echo -e "finalvcf\\tfolder_2\\t变异基因型文件" >> $output_data/WGS.result.info.txt
#echo -e "pop.snpindel.final.vcf.gz\\tfile\\t所有变异文件(包含snp和indel)" >> $output_data/WGS.result.info.txt
#echo -e "pop.snpindel.final.vcf.gz.tbi\\tfile\\t所有变异文件(包含snp和indel)的索引文件" >> $output_data/WGS.result.info.txt
#echo -e "snp\\tfolder_2\\t变异基因型(snp)文件" >> $output_data/WGS.result.info.txt
#echo -e "snpEFF.snp.genes.txt\\tfile\\t使用snpEFF注释变异基因型(snp)文件的基因总表" >> $output_data/WGS.result.info.txt
#echo -e "snpEFF.snp.csv\\tfile\\t变异基因型(snp)文件注释统计总表" >> $output_data/WGS.result.info.txt
#echo -e "snp.summary.html\\tfile\\t变异基因型(snp)文件总结" >> $output_data/WGS.result.info.txt
#echo -e "snp_type_distribution.txt\\tfile\\t变异基因型(snp)文件snp类型分布统计表" >> $output_data/WGS.result.info.txt
#echo -e "snp.stat\\tfile\\t变异基因型(snp)文件单样品snp统计表" >> $output_data/WGS.result.info.txt
#echo -e "snp.stat.all\\tfile\\t变异基因型(snp)文件所有样品snp统计表总表" >> $output_data/WGS.result.info.txt
#echo -e "snp_anno.stat\\tfile\\t变异基因型(snp)文件注释统计表表" >> $output_data/WGS.result.info.txt
#echo -e "pop_snp_type_distribution.pdf/png\\tfile\\t变异基因型(snp)文件snp类型分布统计图" >> $output_data/WGS.result.info.txt
#echo -e "indel\\tfolder_2\\t变异基因型(indel)文件" >> $output_data/WGS.result.info.txt
#echo -e "snpEFF.indel.genes.txt\\tfile\\t使用indelEFF注释变异基因型(indel)文件的基因总表" >> $output_data/WGS.result.info.txt
#echo -e "snpEFF.indel.csv\\tfile\\t变异基因型(indel)文件注释统计总表" >> $output_data/WGS.result.info.txt
#echo -e "indel.summary.html\\tfile\\t变异基因型(indel)文件总结" >> $output_data/WGS.result.info.txt
#echo -e "indel_type_distribution.txt\\tfile\\t变异基因型(indel)文件indel类型分布统计表" >> $output_data/WGS.result.info.txt
#echo -e "indel.stat\\tfile\\t变异基因型(indel)文件单样品indel统计表" >> $output_data/WGS.result.info.txt
#echo -e "indel_anno.stat\\tfile\\t变异基因型(indel)文件注释统计表表" >> $output_data/WGS.result.info.txt
#echo -e "pop_indel_type_distribution.pdf/png\\tfile\\t变异基因型(indel)文件indel类型分布统计图" >> $output_data/WGS.result.info.txt
#
#############05
#echo -e "05.anno\\tfolder_1\\t基因注释信息结果汇总" >> $output_data/WGS.result.info.txt
#echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#echo -e "pop.summary.xls\\tfile\\t变异注释信息结果样品汇总统计表" >> $output_data/WGS.result.info.txt
#echo -e "pop.stat.xls\\tfile\\t变异注释信息结果各样品统计表" >> $output_data/WGS.result.info.txt
#echo -e "pop.eggnog.stat.xls\\tfile\\t变异注释信息结果eggnog方法统计表" >> $output_data/WGS.result.info.txt
#echo -e "pop.go.stat.xls\\tfile\\t变异注释信息结果go方法统计表" >> $output_data/WGS.result.info.txt
#echo -e "pop.kegg.stat.xls\\tfile\\t变异注释信息结果kegg方法统计表" >> $output_data/WGS.result.info.txt
#echo -e "pop.pfam.stat.xls\\tfile\\t变异注释信息结果pfam方法统计表" >> $output_data/WGS.result.info.txt
#echo -e "kegg.enrichment.xls\\tfile\\t变异检测kegg富集分析结果表" >> $output_data/WGS.result.info.txt
#echo -e "kegg.enrichment.pdf\\tfile\\t变异检测kegg富集分析图" >> $output_data/WGS.result.info.txt
#echo -e "kegg.enrichment.png\\tfile\\t变异检测kegg富集分析图" >> $output_data/WGS.result.info.txt
#echo -e "go.enrichment.xls\\tfile\\t变异检测go富集分析结果表" >> $output_data/WGS.result.info.txt
#echo -e "go.enrichment.pdf\\tfile\\t变异检测go富集分析图" >> $output_data/WGS.result.info.txt
#echo -e "go.enrichment.png\\tfile\\t变异检测go富集分析图" >> $output_data/WGS.result.info.txt
#
###############06
#echo -e "06.circos\\tfolder_1\\t变异检测结果可视化圈图" >> $output_data/WGS.result.info.txt
#echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#echo -e "circos.png\\tfile\\t结果可视化圈图" >> $output_data/WGS.result.info.txt
#echo -e "circos.svg\\tfile\\t结果可视化圈图" >> $output_data/WGS.result.info.txt
#
###############07
#echo -e "07.ssr\\tfolder_1\\tssr统计结果汇总" >> $output_data/WGS.result.info.txt
#echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#echo -e "table_stat.out\\tfile\\tssr统计结果汇总表" >> $output_data/WGS.result.info.txt

##############08
echo -e "03.tdna\\tfolder_1\\t插入位点结果汇总" >> $output_data/WGS.result.info.txt
echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
echo -e "tdna_coverage\\tfolder_2\\t样品插入位点覆盖统计结果汇总" >> $output_data/WGS.result.info.txt
for j in `ls 03.tdna/tdna_coverage/*.png |cut -d "." -f1 |sort -u`;do echo $j file $j样品插入位点覆盖统计图; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
echo -e "aimhii\\tfolder_2\\tAIM-HII插入位点结果汇总" >> $output_data/WGS.result.info.txt
for j in `ls 03.tdna/aimhii/*.csv |cut -d "/" -f3 |sort -u`;do echo $j file $j样本aimhii插入位点结果具体文件; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
echo -e "aimhii.stat.xls\\tfile\\tAIM-HII结果统计表" >> $output_data/WGS.result.info.txt
echo -e "tdnascan\\tfolder_2\\tTDNAscan插入位点结果汇总" >> $output_data/WGS.result.info.txt
echo -e "tdna.stat.xls\\tfile\\tTDNAscan结果统计表" >> $output_data/WGS.result.info.txt

###############SV
##文件夹名字
#if [ -d *.sv ];then
#    echo "sv yes"
#    sv_dir=`find -name *.sv |cut -d "/" -f2`
#    echo -e "$sv_dir\\tfolder_1\\tSV变异结果汇总" >> $output_data/WGS.result.info.txt
#    echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#    echo -e "stat\\tfolder_2\\tSV变异结果统计汇总" >> $output_data/WGS.result.info.txt
#    echo -e "stat.txt\\tfile\\t所有样本SV变异结果统计汇总表" >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.sv.length.xls file ${j}样本SV变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.BND.length.pdf file ${j}样本SV-BND变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.BND.length.png file ${j}样本SV-BND变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.DEL.length.pdf file ${j}样本SV-DEL变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.DEL.length.png file ${j}样本SV-DEL变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.DUP.length.pdf file ${j}样本SV-DUP变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.DUP.length.png file ${j}样本SV-DUP变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.INS.length.pdf file ${j}样本SV-INS变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.INS.length.png file ${j}样本SV-INS变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.INV.length.pdf file ${j}样本SV-INV变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    for j in `ls 01.fastq_qc/json |cut -d "." -f1 |sort -u`;do echo $j.INV.length.png file ${j}样本SV-INV变异长度统计表; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#    echo -e "vcf\\tfolder_2\\tSV变异文件" >> $output_data/WGS.result.info.txt
#    echo -e "germline_sv.vcf.gz\\tfile\\t所有样品SV变异文件" >> $output_data/WGS.result.info.txt
#    echo -e "germline_sv.vcf.gz.tbi\\tfile\\t所有样品SV变异文件索引文件" >> $output_data/WGS.result.info.txt
#    echo "sv echo yes"
#fi
#
#
###############CNV
#if [ -d *.cnv ];then
#    cnv_dir=`find -name *.cnv |cut -d "/" -f2`
#    echo -e "$cnv_dir\\tfolder_1\\tCNV变异结果汇总" >> $output_data/WGS.result.info.txt
#    echo -e "结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/WGS.result.info.txt
#    echo -e "all.cnvkit.stats.xls\\tfile\\t所有样本CNV变异长度统计表" >> $output_data/WGS.result.info.txt
#    echo -e "pop.cnv.merge.vcf.gz\\tfile\\t所有样品CNV变异文件" >> $output_data/WGS.result.info.txt
#    echo -e "pop.cnv.merge.vcf.gz.tbi\\tfile\\t所有样品CNV变异文件索引文件" >> $output_data/WGS.result.info.txt
#    for j in `ls $cnv_dir |grep -v "pop" |grep -v "all" |grep -v "结果说明文档" |cut -d "." -f1 |sort -u`;do echo -e "${j} folder_2 ${j}样本CNV变异详情\\ncnvkit.DEL.pdf file ${j}样品CNV-缺失类型的图\\ncnvkit.DEL.png file ${j}样品CNV-缺失类型的图\\ncnvkit.DUP.pdf file ${j}样品CNV-重复类型的图\\ncnvkit.DUP.png file ${j}样品CNV-重复类型的图\\ncnvkit.length.xls file ${j}样品CNV变异长度统计表\\ncnvkit.stat.xls file ${j}样本CNV变异结果统计表"; done |sed 's/ /\t/g' >> $output_data/WGS.result.info.txt
#fi
#
#
#