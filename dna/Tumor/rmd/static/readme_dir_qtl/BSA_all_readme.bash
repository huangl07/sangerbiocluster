echo "01.vcf2table	folder_1	变异检测分析结果
01.vcf2table_readme.txt	file	readme文件
pop.table	file	变异检测统计列表
pop.final.anno.xls	file	变异检测带基因注释表统计列表
snp_indel_gene.stat.xls	file	变异检测统计列表统计
02.ED	folder_1	ED方法获取的定位结果文件夹
02.ED_readme.txt	file	readme文件
pop.ED.sliding.result.xls	file	ED方法结果文件
ED.region_stat.xls	file	ED方法定位区域统计表
ED.all.table.xls	file	ED方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "pop.ED.$i.index.png	file	不同染色体上ED方法的曼哈顿图（png格式）
pop.ED.$i.index.pdf	file	不同染色体上ED方法的曼哈顿图（pdf格式）" >> BSA.info.result
done

echo "pop.ED.index.pdf	file	ED方法的整体曼哈顿图（pdf格式）
pop.ED.index.png	file	ED方法的整体曼哈顿图（png格式）
03.index	folder_1	index方法获取的定位结果文件夹
03.index_readme.txt	file	readme文件
pop.index.bootstrap.result.xls	file	index方法结果文件
index.region_stat.xls	file	index方法定位区域统计表
index.all.table.xls	file	index方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "pop.index.*.index.pdf	file	不同染色体上index方法的曼哈顿图（pdf格式）
pop.index.*.index.png	file	不同染色体上index方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "pop.index.index.pdf	file	index方法的整体曼哈顿图（pdf格式）
pop.index.index.png	file	index方法的整体曼哈顿图（png格式）
04.Gprime	folder_1	Gprime方法获取的定位结果文件夹
04.Gprime_readme.txt	file	readme文件
pop.Gprime.result.xls	file	Gprime方法统计结果文件
Gprime.region_stat.xls	file	Gprime方法定位区域位点和基因详情表
Gprime.all.table.xls	file	Gprime方法的整体曼哈顿图（png格式）" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "Gprime.*.index.pdf	file	不同染色体上Gprime方法的曼哈顿图（pdf格式）
Gprime.*.index.png	file	不同染色体上Gprime方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "Gprime.index.pdf	file	Gprime方法的整体曼哈顿图（pdf格式）
Gprime.index.png	file	Gprime方法的整体曼哈顿图（png格式）
05.loess	folder_1	loess方法获取的定位结果文件夹
05.loess_readme.txt	file	readme文件
pop.loess.bootstrap.result.xls	file	loess方法统计结果文件
loess.region_stat.xls	file	loess方法定位区域统计表
loess.all.table.xls	file	loess方法定位区域位点和基因详情表" >> BSA.info.result

for i in $(cut -f1 $workdir/data_info/chr.list);do
echo "loess.*.index.pdf	file	不同染色体上loess方法的曼哈顿图（pdf格式）
loess.*.index.pdf	file	不同染色体上loess方法的曼哈顿图（png格式）" >> BSA.info.result
done

echo "loess.index.pdf	file	loess方法的整体曼哈顿图（pdf格式）
loess.index.png	file	loess方法的整体曼哈顿图（png格式）
06.enrich	folder_1	候选基因GO和KEGG富集分析
06.enrich_readme.txt	file	readme文件
index	folder_2	候选基因GO和KEGG富集分析（index方法）
index.genes_abstract.list	file	定位区域基因bed文件
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
ED	folder_2	候选基因GO和KEGG富集分析（ED方法）
ED.genes_abstract.list	file	定位区域基因bed文件
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
Gprime.genes_abstract.list	file	定位区域基因bed文件
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
loess	folder_2	候选基因GO和KEGG富集分析（loess方法）
loess.genes_abstract.list	file	定位区域基因bed文件
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
" >> BSA.info.result
