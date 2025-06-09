#GMAP路径 $1 如：/mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/
#output路径 $2
#data_tmp结果路径 $3 data_tmp结果路径 如/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/output/tmp/
#GMAP结果路径 $5
Usage(){
    echo "Usage:"
    echo "bash run_report.sh -m gmap_path -g gmap_result -o output_path -w data_tmp -s wgs_result"
    echo "-m GMAP路径 如：/mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/"
    echo "-g GMAP结果路径"
    echo "-o output路径"
    echo "-d data/下面保存project.info/info.log"
    exit -1
}
while getopts ":g:o:w:m:" opt_name
do
    case $opt_name in
        m) gmap_path="$OPTARG"
            ;;
        g) gmap_result="$OPTARG"
            ;;
        o) output_path="$OPTARG"
            ;;
        d) data_tmp="$OPTARG"
            ;;
        :) Usage;;
    esac
done

source ~/app/bioinfo/dna/new.rc
cd ${output_path}
#Rscript ${gmap_path}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path}/遗传图谱.pdf --rmd ${gmap_path}/report_tmpl/report.rmd --wgs_result ${data_tmp} --gmap_result ${gmap_result}
#Rscript ${gmap_path}/report_tmpl/rmarkdown.R --format html --data_release ${output_path}/遗传图谱.html --rmd ${gmap_path}/report_tmpl/report.rmd --wgs_result ${data_tmp} --gmap_result ${gmap_result}
Rscript ${gmap_path}/report_tmpl/rmarkdown.R --format pdf  --outfile ${output_path}/遗传图谱.pdf --data_release ${output_path} --rmd ${gmap_path}/report_tmpl/report.rmd --data_result ${data_tmp} --gmap_result ${gmap_result}
Rscript ${gmap_path}/report_tmpl/rmarkdown.R --format html --outfile ${output_path}/遗传图谱.html --data_release ${output_path} --rmd ${gmap_path}/report_tmpl/report.rmd --data_result ${data_tmp} --gmap_result ${gmap_result}
#示例代码
#Rscript /mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/report_tmpl/rmarkdown.R --format pdf --outfile /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/report/遗传图谱.html --rmd  /mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/report_tmpl/report.rmd --wgs_result /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/output/tmp/ --gmap_result /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/ --data_release /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/report_new --stat_path /mnt/lustre/sanger_workspaceWgsV4/20231222/WgsV4_a2o2_0flc2pd5np50mqc18nuipu
#Rscript /mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/report_tmpl/rmarkdown.R --format pdf --outfile /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/report/遗传图谱.pdf --rmd  /mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/report_tmpl/report.rmd --wgs_result /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/output/tmp/ --gmap_result /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/ --data_release /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/report_new --stat_path /mnt/lustre/sanger_workspaceWgsV4/20231222/WgsV4_a2o2_0flc2pd5np50mqc18nuipu
