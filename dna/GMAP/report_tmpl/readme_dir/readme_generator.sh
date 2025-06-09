output_data=$1
touch $output_data/GMAP.info.txt
echo -e "01.vcf_filter\\tfolder_1\\t遗传标记筛选及分析结果\\n01.vcf_filter_结果说明文档.txt\\tfile\\treadme文件\\npop.snp_anno.xls\\tfile\\tSNP功效信息统计表\\npop.indel_anno.xls\\tfile\\tINDEL功效信息统计表\\ntotal.filtered.stat.xls\\tfile\\tINDEL功效信息统计表\\n02.genetic_map\\tfolder_1\\t遗传图谱构建结果\\n02.genetic_map_结果说明文档.txt\\tfile\\t文件详细介绍说明\\ntotal.mapstat.xls\\tfile\\tSNP功效信息统计表\\ntotal.map.pdf\\tfile\\t遗传图谱示意图(pdf版本)\\ntotal.map.png\\tfile\\t遗传图谱示意图(png版本)\\n03.evaluate\\tfolder_1\\t图谱评估结果\\n03.evaluate_结果说明文档.txt\\tfile\\t文件详细介绍说明\\ntotal.bin.pdf\\tfile\\t连锁群的单体来源图(pdf版本)\\ntotal.bin.png\\tfile\\t连锁群的单体来源图(png版本)\\ntotal.phy.pdf\\tfile\\t遗传图谱与基因组线性图(pdf版本)
total.phy.png\\tfile\\t遗传图谱与基因组线性图(pdf版本)\\ntotal.phy.spearman.xls\\tfile\\t共线性统计表" >> $output_data/GMAP.info.txt
# 判断是否存在 *.sexAver 文件
if ls 03.evaluate/*heatMap.sexAver.* 1> /dev/null 2>&1; then
    for j in `ls 03.evaluate/*heatMap.sexAver.png | cut -d "/" -f2 | cut -d "." -f1 | sort -u`; do 
        echo -e "${j}.heatMap.sexAver.png\\tfile\\t${j}重组关系热图，图片上方的编号代表染色体号(png版本)" >> $output_data/GMAP.info.txt
    done
    
    for j in `ls 03.evaluate/*heatMap.sexAver.pdf | cut -d "/" -f2 |  cut -d "." -f1 |sort -u`; do
        echo -e "${j}.heatMap.sexAver.pdf\\tfile\\t${j}重组关系热图，图片上方的编号代表染色体号(pdf版本)" >> $output_data/GMAP.info.txt
    done
fi

# 判断是否存在 *.heatMap 文件
if ls 03.evaluate/*heatMap.* 1> /dev/null 2>&1; then
    for j in `ls 03.evaluate/*heatMap.png | cut -d "/" -f2 | cut -d "." -f1 | sort -u`; do 
        echo -e "${j}.heatMap.png\\tfile\\t${j}重组关系热图，图片上方的编号代表染色体号(png版本)" >> $output_data/GMAP.info.txt
    done
    
    for j in `ls 03.evaluate/*heatMap.pdf | cut -d "/" -f2 |  cut -d "." -f1 |sort -u`; do
        echo -e "${j}.heatMap.pdf\\tfile\\t${j}重组关系热图，图片上方的编号代表染色体号(pdf版本)" >> $output_data/GMAP.info.txt
    done
fi

if [ -d 04.qtl ];then
    echo -e "04.qtl\\tfolder_1\\tqtl结果目录" >> $output_data/GMAP.info.txt
    echo -e "04.qtl_结果说明文档.txt\\tfile\\t文件详细介绍说明" >> $output_data/GMAP.info.txt
    echo -e "drug.stat.xls\\tfile\\t药物靶点统计表" >> $output_data/GMAP.info.txt
    for j in `ls 04.qtl| grep -v "结果说明文档" |sort -u`;do 
        echo -e "${j}\\tfolder_2\\t${j}qtl结果文件" >> $output_data/GMAP.info.txt
        echo -e "enrich\\tfolder_3\\t${j}qtl性状的enrich结果文件">> $output_data/GMAP.info.txt
        echo -e "${j}.region.gene.xls\\tfile\\t定位区域内基因统计表" >> $output_data/GMAP.info.txt
        echo -e "${j}.region.variant.xls\\tfile\\t定位区域内突变位点统计表" >> $output_data/GMAP.info.txt
        echo -e "GO_result\\tfolder_4\\t高频突变基因结果目录" >> $output_data/GMAP.info.txt
        echo -e "${j}_GOenrichment.xls\\tfile\\tGO富集结果表(all)" >> $output_data/GMAP.info.txt
        echo -e "${j}_GOenrichment_0.05.xls\\tfile\\tGO富集结果表(显著阈值在0.05之下)" >> $output_data/GMAP.info.txt
        echo -e "${j}_GOenrichment.pdf\\tfile\\tGO富集结果图" >> $output_data/GMAP.info.txt
        echo -e "${j}_GOenrichment.png\\tfile\\tGO富集结果图" >> $output_data/GMAP.info.txt
        echo -e "KEGG_result\\tfolder_4\\t高频突变基因结果目录" >> $output_data/GMAP.info.txt
        echo -e "${j}_KEGGenrichment.xls\\tfile\\tKEGG富集结果表(all)" >> $output_data/GMAP.info.txt
        echo -e "${j}_KEGGenrichment_0.05.xls\\tfile\\tKEGG富集结果表(显著阈值在0.05之下)" >> $output_data/GMAP.info.txt
        echo -e "${j}_KEGGenrichment.pdf\\tfile\\tKEGG富集结果图" >> $output_data/GMAP.info.txt
        echo -e "${j}_KEGGenrichment.png\\tfile\\tKEGG富集结果图" >> $output_data/GMAP.info.txt
        echo -e "qtl\\tfolder_3\\t${j}qtl性状结果文件" >> $output_data/GMAP.info.txt
        echo -e "${j}.qtl-result.result\\tfile\\t性状定位结果表" >> $output_data/GMAP.info.txt
        echo -e "${j}.detail.result\\tfile\\t性状定位结果详细结果表" >> $output_data/GMAP.info.txt
        echo -e "${j}.pheno.png/pdf\\tfile\\t性状频率分布图(png/pdf版本)" >> $output_data/GMAP.info.txt
        echo -e "${j}.png/pdf\\tfile\\tQTL定位曼哈顿图(png/pdf版本)" >> $output_data/GMAP.info.txt
        echo -e "${j}.scan.png/pdf\\tfile\\tQTL定位结果图(png/pdf版本)" >> $output_data/GMAP.info.txt
        echo -e "${j}.pxg.pdf\\tfile\\t区域内不同基因型的性状分布图(pdf版本,如果定位结果不存在则没有)" >> $output_data/GMAP.info.txt
        echo -e "${j}.pxg.png\\tfile\\t区域内不同基因型的性状分布图(png版本,如果定位结果不存在则没有)" >> $output_data/GMAP.info.txt
    done
fi
