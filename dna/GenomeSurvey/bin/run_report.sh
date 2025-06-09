#! bash -ue

if [ $# -ne 4 ]; then
    echo 'Usage: '$(basename $0)' gs_result output_dir report_template project.info'
    echo 'gs_result: genome_survey输出目录'
    echo 'output_dir: 最终目录'
    echo 'report_template: 报告模板'
    echo 'project.info: 项目信息'
    exit 1
fi

rm -rf report
mkdir -p report
mkdir -p $2/data
cp -rl $1/qc $2/data/01.fastq_qc
cp -rl $1/gs $2/data/02.genome_survey
rm -f report/gs_result && ln -s $(readlink -f $2)/data report/wgs_result
cp $4 report/project.info
cd report
cp -rf $3/template/* .
Rscript rmarkdown.R --rmd report.rmd --format html --outfile 基因组调查分析报告.html --type genome_survey
Rscript rmarkdown.R --rmd report.rmd --format pdf --outfile 基因组调查分析报告.pdf --type genome_survey
cd -
cp -l report/基因组调查分析报告.{html,pdf} $2
