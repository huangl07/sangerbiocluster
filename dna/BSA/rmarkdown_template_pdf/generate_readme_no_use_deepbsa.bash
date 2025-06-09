echo "01.vcf2table	folder_1	变异检测分析结果
01.vcf2table_结果说明文档.txt	file	readme文件
pop.table	file	变异检测统计列表
pop.final.anno.xls	file	变异检测带基因注释表统计列表
snp_indel_gene.stat.xls	file	变异检测统计表
snp_anno.xls	file	SNP功效信息统计表
indel_anno.xls	file	InDel功效信息统计表
02.index	folder_1	index方法获取的定位结果文件夹
02.index_readme.txt	file	readme文件
pop.index.bootstrap.result.xls	file	index方法结果文件
index.region_stat.xls	file	index方法定位区域统计表
index.all.table.xls	file	index方法定位区域位点和基因详情表
index.detail.result.xls	file	index计算详细结果表" > $1/BSA.info.result

for i in $(cut -f1 $2);do
echo "index.$i.png	file	不同染色体上index方法的曼哈顿图（png格式）
index.$i.pdf	file	不同染色体上index方法的曼哈顿图（pdf格式）" >> $1/BSA.info.result
done

echo "index.png	file	index方法的整体曼哈顿图（png格式）
index.pdf	file	index方法的整体曼哈顿图（pdf格式）
03.loess	folder_1	loess方法获取的定位结果文件夹
03.loess_结果说明文档.txt	file	readme文件
pop.loess.bootstrap.result.xls	file	loess方法结果文件
loess.region_stat.xls	file	loess方法定位区域统计表
loess.all.table.xls	file	loess方法定位区域位点和基因详情表
loess.detail.result.xls	file	loess计算详细结果表" >> $1/BSA.info.result

for i in $(cut -f1 $2);do
echo "loess.$i.pdf	file	不同染色体上loess方法的曼哈顿图（pdf格式）
loess.$i.png	file	不同染色体上loess方法的曼哈顿图（png格式）" >> $1/BSA.info.result
done

echo "loess.pdf	file	loess方法的整体曼哈顿图（pdf格式）
loess.png	file	loess方法的整体曼哈顿图（png格式）
04.ED	folder_1	ED方法获取的定位结果文件夹
04.ED_结果说明文档.txt	file	readme文件
pop.ED.sliding.result.xls	file	ED方法统计结果文件
ED.region_stat.xls	file	ED方法定位区域位点和基因详情表
ED.all.table.xls	file	ED方法的整体曼哈顿图（png格式）
ED.detail.result.xls	file	ED计算详细结果表" >> $1/BSA.info.result

for i in $(cut -f1 $2);do
echo "ED.$i.pdf	file	不同染色体上ED方法的曼哈顿图（pdf格式）
ED.$i.png	file	不同染色体上ED方法的曼哈顿图（png格式）" >> $1/BSA.info.result
done

echo "ED.pdf	file	ED方法的整体曼哈顿图（pdf格式）
ED.png	file	ED方法的整体曼哈顿图（png格式）
05.Gprime	folder_1	Gprime方法获取的定位结果文件夹
05.Gprime_结果说明文档.txt	file	readme文件
pop.Gprime.result.xls	file	Gprime方法统计结果文件
Gprime.region_stat.xls	file	Gprime方法定位区域统计表
Gprime.all.table.xls	file	Gprime方法定位区域位点和基因详情表
Gprime.detail.result.xls	file	Gprime计算详细结果表" >> $1/BSA.info.result

for i in $(cut -f1 $2);do
echo "Gprime.$i.pdf	file	不同染色体上Gprime方法的曼哈顿图（pdf格式）
Gprime.$i.pdf	file	不同染色体上Gprime方法的曼哈顿图（png格式）" >> $1/BSA.info.result
done

echo "Gprime.pdf	file	Gprime方法的整体曼哈顿图（pdf格式）
Gprime.png	file	Gprime方法的整体曼哈顿图（png格式）
06.enrich	folder_1	候选基因GO和KEGG富集分析
06.enrich_结果说明文档.txt	file	readme文件
region	folder_2	各方法定位区域注释基因信息表
index.anno.xls	file	index定位区域注释基因信息表
Gprime.anno.xls	file	Gprime定位区域注释基因信息表
ED.anno.xls	file	ED定位区域注释基因信息表
loess.anno.xls	file	loess定位区域注释基因信息表
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
07.primer_design	folder_1	候选区域引物设计结果
07.primer_design_结果说明文档.txt	file	readme文件
index	folder_2	候选区域引物设计结果（index方法）
index.sanger.snp.result.xls	file	常规变异位点引物设计(snp)
index.sanger.indel.result.xls	file	常规变异位点引物设计(indel)
index.caps.xls	file	CAPS引物设计结果
index.dcaps.xls	file	dCAPS引物设计结果
index.kasp.result.xls file	Kasp引物设计结果
loess	folder_2	候选区域引物设计结果（loess方法）
loess.sanger.snp.result.xls	file	常规变异位点引物设计(snp)
loess.sanger.indel.result.xls	file	常规变异位点引物设计(indel)
loess.caps.result.xls	file	CAPS引物设计结果
loess.dcaps.result.xls	file	dCAPS引物设计结果
loess.kasp.result.xls	file	Kasp引物设计结果
ED	folder_2	候选区域引物设计结果（ED方法）
ED.sanger.snp.result.xls	file	常规变异位点引物设计(snp)
ED.sanger.indel.result.xls	file	常规变异位点引物设计(indel)
ED.caps.result.xls	file	CAPS引物设计结果
ED.dcaps.result.xls	file	dCAPS引物设计结果
ED.kasp.result.xls	file	Kasp引物设计结果
Gprime	folder_2	候选区域引物设计结果（Gprime方法）
Gprime.sanger.snp.result.xls	file	常规变异位点引物设计(snp)
Gprime.sanger.snp.result.xls	file	常规变异位点引物设计(indel)
Gprime.caps.result.xls	file	CAPS引物设计结果
Gprime.dcaps.result.xls	file	dCAPS引物设计结果
Gprime.kasp.result.xls	file	Kasp引物设计结果" >> $1/BSA.info.result
