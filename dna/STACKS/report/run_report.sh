Usage(){
    echo "Usage:"
    echo "bash run_report.sh -r result_dir -o output_realease_dir -w wgs_result_dir -s wgs_result -m report_path"
    echo "-r 无参结果路径"
    echo "-o output路径"
    echo "-m report_path路径 如 /mnt/ilustre/users/dandan.zhang/05_project_2024/070_MJ20240515344_FX2024071800414_tags/dna/STACKS"
    echo "-p project_info文件 请自己根据wgs生成的project.info生成类似的，示例文件如下：/mnt/ilustre/users/dandan.zhang/05_project_2024/070_MJ20240515344_FX2024071800414_tags/project.info"
    exit -1
}
while getopts ":r:o:p:m:" opt_name
do
    case $opt_name in
        r) result_dir="$OPTARG"
            ;;
        o) realease_dir="$OPTARG"
            ;;
        p) project_info="$OPTARG"
            ;;
        m) rmd_dir="$OPTARG"
            ;;
        :) Usage;;
    esac
done
mkdir -p ${realease_dir}
mkdir -p ${realease_dir}/info
mkdir -p ${realease_dir}/data_release
mkdir -p ${realease_dir}/data_release/01.stat
mkdir -p ${realease_dir}/data_release/02.variants


cp -r ${rmd_dir}/report ${realease_dir}
mkdir -p ${realease_dir}/report/info
info_dir=`readlink -f ${realease_dir}/report/info`

grep -v "rawdata" ${result_dir}/09.results/pop.qc | cut -f 1-12 > ${realease_dir}/qc.temp
cat ${rmd_dir}/report/qc_header.txt ${realease_dir}/qc.temp > ${realease_dir}/data_release/01.stat/qc.stat.xls
rm -rf ${realease_dir}/qc.temp
cp ${result_dir}/09.results/sample.stat ${realease_dir}/data_release/01.stat/sample_tags.stat.xls
cp ${result_dir}/09.results/populations.loci.fa ${realease_dir}/data_release/02.variants/ref.loci.fa
cp ${result_dir}/09.results/populations.loci.fa.fai ${realease_dir}/data_release/02.variants/ref.loci.fa.fai
cp ${result_dir}/09.results/populations.final.vcf.gz ${realease_dir}/data_release/02.variants/populations.final.vcf.gz
cp ${result_dir}/09.results/populations.final.vcf.gz.tbi ${realease_dir}/data_release/02.variants/populations.final.vcf.gz.tbi
cp -r ${result_dir}/09.result/fig ${realease_dir}/data_release/01.stat/
cp ${result_dir}/09.results/snp.stat.xls ${realease_dir}/data_release/02.variants/
cp ${result_dir}/09.results/snp.stat.all.xls ${realease_dir}/data_release/02.variants/

cp ${rmd_dir}/report/readme/01.stat.readme.txt ${realease_dir}/data_release/01.stat/结果说明文档.txt
cp ${rmd_dir}/report/readme/02.variants.readme.txt ${realease_dir}/data_release/02.variants/结果说明文档.txt

cp ${project_info} ${info_dir}/project.info
cp ${realease_dir}/data_release/01.stat/qc.stat.xls ${info_dir}/qc.stat
cp ${result_dir}/09.results/sample.stat ${info_dir}/sample.stat
cp ${result_dir}/09.result/fig/*png ${info_dir}
cp ${result_dir}/09.results/snp.stat.xls ${info_dir}

echo "Rscript ${realease_dir}/report/rmarkdown.R --rmd ${realease_dir}/report/report_noref.rmd --format html --outfile ${realease_dir}/简化基因组（无参）报告.html"
echo "Rscript ${realease_dir}/report/rmarkdown.R --format pdf --outfile ${realease_dir}/简化基因组（无参）报告.pdf --rmd ${realease_dir}/report/report_noref.rmd"
Rscript ${realease_dir}/report/rmarkdown.R --rmd ${realease_dir}/report/report_noref.rmd --format html --outfile ${realease_dir}/简化基因组（无参）报告.html
Rscript ${realease_dir}/report/rmarkdown.R --format pdf --outfile ${realease_dir}/简化基因组（无参）报告.pdf --rmd ${realease_dir}/report/report_noref.rmd

bash ${realease_dir}/report/readme/readme_generator.sh ${realease_dir}/data_release

cd ${realease_dir}/data_release/
python3 ${realease_dir}/report/readme/readme.py -f ${realease_dir}/data_release/noref.result.info.txt -t RAD_noref -n RAD_noref -b ${realease_dir}/report/readme/
