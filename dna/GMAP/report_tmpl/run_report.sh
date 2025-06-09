#GMAP路径 $1 如：/mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/
#output路径 $2
#wgs_tmp结果路径 $3 wgs_tmp结果路径 如/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/output/tmp/
#GMAP结果路径 $5
Usage(){
    echo "Usage:"
    echo "bash run_report.sh -m gmap_path -g gmap_result -o output_path -w wgs_result"
    echo "-m GMAP路径 如：/mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/"
    echo "-g GMAP结果路径"
    echo "-o output路径"
    echo "-w wgs_result结果路径 如/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/"
    echo "-t type,RAD/RAD_noref/WGS,请选择模式"
    exit -1
}
while getopts ":g:o:w:m:t:" opt_name
do
    case $opt_name in
        m) gmap_path="$OPTARG"
            ;;
        g) gmap_result="$OPTARG"
            ;;
        o) output_path="$OPTARG"
            ;;
        w) wgs_result="$OPTARG"
            ;;
        t) type="$OPTARG"
            ;;
        :) Usage;;
    esac
done
source ~/app/bioinfo/dna/new.rc
mkdir -p ${output_path}/data_release
Result_dir=${output_path}/data_release
cd ${Result_dir}
pwd

mkdir -p $Result_dir/01.vcf_filter
mkdir -p $Result_dir/02.genetic_map
mkdir -p $Result_dir/03.evaluate


if [[ "$type" == "WGS" ]]; then
    cp -r ${gmap_path}/report_tmpl $Result_dir
    cp ${gmap_result}/08.evaluate/total.filtered.stat $Result_dir/01.vcf_filter/total.filtered.stat.xls
    cp ${wgs_result}/output/published/data/04.snpIndel/indel/indel_anno.stat.xls $Result_dir/01.vcf_filter/pop.snp_anno.xls
    cp ${wgs_result}/output/published/data/04.snpIndel/snp/snp_anno.stat.xls $Result_dir/01.vcf_filter/pop.indel_anno.xls
    cp ${gmap_result}/08.evaluate/total.map.png $Result_dir/02.genetic_map/total.map.png
    cp ${gmap_result}/08.evaluate/total.map.pdf $Result_dir/02.genetic_map/total.map.pdf
    cp ${gmap_result}/08.evaluate/total.mapstat $Result_dir/02.genetic_map/total.mapstat.xls

    cp ${gmap_result}/09.result/*.heatMap*.pdf $Result_dir/03.evaluate/
    cp ${gmap_result}/09.result/*.heatMap*.png $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.pdf $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.png $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.spearman.xls $Result_dir/03.evaluate/

    cp ${Result_dir}/report_tmpl/readme_dir/01.vcf_filter_结果说明文档.txt $Result_dir/01.vcf_filter/
    cp ${Result_dir}/report_tmpl/readme_dir/02.genetic_map_结果说明文档.txt $Result_dir/02.genetic_map
    cp ${Result_dir}/report_tmpl/readme_dir/03.evalute_结果说明文档.txt $Result_dir/03.evaluate/

    # 检查 qtl_result 文件夹是否存在
    if [ -d "${gmap_result}/qtl_result" ]; then
        run_qtl="yes"
        # 如果存在，则创建 04.qtl 目录（如果它还不存在）
        mkdir -p "$Result_dir/04.qtl"
        # 将 qtl_result/result 目录复制到 $Result_dir/04.qtl/
        cp -rl ${gmap_result}/qtl_result/result/* $Result_dir/04.qtl/
        cp ${Result_dir}/report_tmpl/readme_dir/04.qtl_结果说明文档.txt $Result_dir/04.qtl
        bash ${Result_dir}/report_tmpl/readme_dir/readme_generator.sh ./
        python3 ${Result_dir}/report_tmpl/readme_dir/ReadMe_result.py -f GMAP.info.txt -t GMAP_qtl -n GMAP
    else
        echo "qtl_result 文件夹不存在！运行无qtl步骤"
        run_qtl="no"
        bash ${Result_dir}/report_tmpl/readme_dir/readme_generator.sh ./
        python3 ${Result_dir}/report_tmpl/readme_dir/ReadMe_result.py -f GMAP.info.txt -t GMAP -n GMAP
    fi

    echo "Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type wgs"
    Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type wgs
    echo "Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format html --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type wgs"
    Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format html --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type wgs

    mv ${Result_dir}/report_tmpl/report.html ${Result_dir}/遗传图谱报告.html
    mv ${Result_dir}/report_tmpl/report.pdf ${Result_dir}/遗传图谱报告.pdf
elif [[ "$type" == "RAD" ]]; then
    cp -r ${gmap_path}/report_tmpl $Result_dir
    cp ${gmap_result}/08.evaluate/total.filtered.stat $Result_dir/01.vcf_filter/total.filtered.stat.xls
    cp ${wgs_result}/output/published/data/04.snpIndel/indel/indel_anno.stat.xls $Result_dir/01.vcf_filter/pop.snp_anno.xls
    cp ${wgs_result}/output/published/data/04.snpIndel/snp/snp_anno.stat.xls $Result_dir/01.vcf_filter/pop.indel_anno.xls

    cp ${gmap_result}/08.evaluate/total.map.png $Result_dir/02.genetic_map/total.map.png
    cp ${gmap_result}/08.evaluate/total.map.pdf $Result_dir/02.genetic_map/total.map.pdf
    cp ${gmap_result}/08.evaluate/total.mapstat $Result_dir/02.genetic_map/total.mapstat.xls

    cp ${gmap_result}/08.evaluate/*.heatMap.pdf $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/*.heatMap.png $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.pdf $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.png $Result_dir/03.evaluate/
    cp ${gmap_result}/08.evaluate/total.phy.spearman.xls $Result_dir/03.evaluate/

    cp ${Result_dir}/report_tmpl/readme_dir/01.vcf_filter_结果说明文档.txt $Result_dir/01.vcf_filter/
    cp ${Result_dir}/report_tmpl/readme_dir/02.genetic_map_结果说明文档.txt $Result_dir/02.genetic_map
    cp ${Result_dir}/report_tmpl/readme_dir/03.evalute_结果说明文档.txt $Result_dir/03.evaluate/

    # 检查 qtl_result 文件夹是否存在
    if [ -d "${gmap_result}/qtl_result" ]; then
        run_qtl="yes"
        # 如果存在，则创建 04.qtl 目录（如果它还不存在）
        mkdir -p "$Result_dir/04.qtl"
        # 将 qtl_result/result 目录复制到 $Result_dir/04.qtl/
        cp -rl ${gmap_result}/qtl_result/result/* $Result_dir/04.qtl/
        cp ${Result_dir}/report_tmpl/readme_dir/04.qtl_结果说明文档.txt $Result_dir/04.qtl
        bash ${Result_dir}/report_tmpl/readme_dir/readme_generator.sh ./
        python3 ${Result_dir}/report_tmpl/readme_dir/ReadMe_result.py -f GMAP.info.txt -t GMAP_qtl -n GMAP
    else
        echo "qtl_result 文件夹不存在！运行无qtl步骤"
        run_qtl="no"
        bash ${Result_dir}/report_tmpl/readme_dir/readme_generator.sh ./
        python3 ${Result_dir}/report_tmpl/readme_dir/ReadMe_result.py -f GMAP.info.txt -t GMAP -n GMAP
    fi

    echo "Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type rad"
    Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type rad
    echo "Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format html --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type rad"
    Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format html --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type rad

    mv ${Result_dir}/report_tmpl/report.html ${Result_dir}/遗传图谱报告.html
    mv ${Result_dir}/report_tmpl/report.pdf ${Result_dir}/遗传图谱报告.pdf
elif [[ "$type" == "RAD_noref" ]]; then
    echo "wait"
else
    echo "wait"
fi