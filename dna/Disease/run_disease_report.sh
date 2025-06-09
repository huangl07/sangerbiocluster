#! bash -ue


if [ $# -ne 7 ];
then
    echo 'Usage: '$(basename $0)' wgs_result disease_result output_dir library_type report_template cog asso'
    echo 'wgs_result: wgs结果目录'
    echo 'disease_result: disease输出目录'
    echo 'output_dir: 最终目录'
    echo 'library_type: wgs/wes'
    echo 'report_template: 报告模板'
    echo 'cog: 是否运行共分离'
    echo 'asso: 是否运行疾病基因关联'
    exit 1
fi

WGSRESULT=$(readlink -f $1)
DISRESULT=$(readlink -f $2)
OUTPUTDIR=$(readlink -f $3)

mkdir -p report/tmp
cp -rf ${WGSRESULT}/published/data/* report/tmp
rm -rf report/tmp/05.disease
rm -rf report/tmp/*.html
cp -rl ${DISRESULT}/annotation report/tmp/05.disease
cd report
# python3 $(dirname $0)/bin/file_arrange.py -c $(dirname $0)/bin/disease_file.yml -i tmp -o wgs_result -t
# cat $(dirname $0)/bin/readme/header.html body.html $(dirname $0)/bin/readme/tail.html > 结果目录索引及说明.html
# cp -f 结果目录索引及说明.html $3/data
rm -f wgs_result && ln -s tmp wgs_result
cp -f ${WGSRESULT}/tmp/02.reference/info.log .
cp -f ${WGSRESULT}/tmp/02.reference/project.info .
cp -rf $5/template/* .
rm -f disease_result && ln -s wgs_result disease_result
Rscript rmarkdown.R --rmd report.rmd --format html --outfile 遗传病分析报告.html --params cog=${6} asso=${7} --type ${4}_disease
Rscript rmarkdown.R --rmd report.rmd --format pdf --outfile 遗传病分析报告.pdf --params cog=${6} asso=${7} --type ${4}_disease
rm -rf ${OUTPUTDIR}/data/
mkdir -p ${OUTPUTDIR}/data/
cp -rl wgs_result/* ${OUTPUTDIR}/data/
cp -l 遗传病分析报告.{html,pdf} ${OUTPUTDIR}
