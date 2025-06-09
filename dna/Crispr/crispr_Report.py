# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/06/07
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
from pathlib import Path
import shutil
import os
from bin.bsa_function import run_cmd

def arrange_wgs_v4(sh1):
    with open(sh1, "w") as w:
        ## qc结果整理
        w.write(f"""for i in `cut -f 1 {group_list}`
                    do
                    cp {workflow_result}/published/data/01.fastq_qc/fig/$i.raw.base.png {file_temp}
                    cp {workflow_result}/published/data/01.fastq_qc/fig/$i.clean.base.png {file_temp}
                    cp {workflow_result}/published/data/01.fastq_qc/fig/$i.raw.qual.png {file_temp}
                    cp {workflow_result}/published/data/01.fastq_qc/fig/$i.clean.qual.png {file_temp}
                    done
                """)
        w.write("\n")
        w.write(f"python3 {rmarkdown_scripts_bin}/process_qc_result.py --qc_stat {workflow_result}/published/data/01.fastq_qc/qc.stat.txt \
                --group_file {group_list} --raw_data {file_temp}/rawdata.xls --clean_data {file_temp}/cleandata.xls")
        w.write("\n")
        ## 比对结果统计
        w.write(f"""for i in `cut -f 1 {group_list}`
                    do
                    cp {workflow_result}/published/data/03.mappingStat/$i.genome.coverage.png {file_temp}
                    cp {workflow_result}/published/data/03.mappingStat/$i.insert.png {file_temp}
                    cp {workflow_result}/published/data/03.mappingStat/$i.depth.png {file_temp}
                    done
                """)
        w.write("\n")
        w.write(f"python3 {rmarkdown_scripts_bin}/process_mapping_result.py --all_summary_stats {workflow_result}/published/data/03.mappingStat/all.summary.stats.xls \
                --group_file {group_list} --bsa_align_stat {file_temp}/align_stat.xls")
        w.write("\n")
        ## 变异检测结果整理
        w.write(f"python3 {rmarkdown_scripts_bin}/process_snp_result.py --group_file {group_list} \
                --snp_stat {workflow_result}/published/data/04.snpIndel/snp/snp.stat.xls --snp_stat_bsa {file_temp}/snp_stat.xls")
        w.write("\n")
        w.write(f"python3 {rmarkdown_scripts_bin}/process_indel_result.py --group_file {group_list} \
                --indel_stat {workflow_result}/published/data/04.snpIndel/indel/indel.stat.xls --indel_stat_bsa {file_temp}/indel_stat.xls")
        w.write("\n")
        w.write(f"python3 {rmarkdown_scripts_bin}/process_snp_anno_result.py --group_file {group_list} \
                --snp_anno_stat {workflow_result}/published/data/04.snpIndel/snp/snp_anno.stat.xls \
                --snp_anno_stat_bsa {file_temp}/snp_anno.xls")
        w.write("\n")
        w.write('''cut -f 1,10,12,17-19,21 %s/snp_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
                | sed 's/[ \\t]*$//g' | sed 's/ /\\t/g' > %s/snp_anno_T.xls''' % (file_temp, file_temp))
        w.write("\n")
        w.write(f"cut -f 1,5-8 {file_temp}/snp_anno.xls > {file_temp}/snp_impact.xls")
        w.write("\n")
        w.write(f"python3 {rmarkdown_scripts_bin}/process_indel_anno_result.py --group_file {group_list} \
                --indel_anno_stat {workflow_result}/published/data/04.snpIndel/indel/indel_anno.stat.xls \
                --indel_anno_stat_bsa {file_temp}/indel_anno.xls")
        w.write("\n")
        w.write('''cut -f 1,16,17,19,27,29,30 %s/indel_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
                | sed 's/[ \\t]*$//g' | sed 's/ /\\t/g' > %s/indel_anno_T.xls''' % (file_temp, file_temp))
        w.write("\n")
        w.write(f"cut -f 1,6-9 {file_temp}/indel_anno.xls > {file_temp}/indel_impact.xls")
        w.write("\n")
        w.write(f"cp {workflow_result}/tmp/02.reference/project.info {file_temp}")
        w.write("\n")
        w.write(f"cp {workflow_result}/tmp/02.reference/info.log {file_temp}")
        w.write("\n")
        w.write(f"cp {group_list} {file_temp}/group_list.txt")

def arrange_crispr_off_target(sh2):
    with open(sh2, "w") as w:
        w.write(f"cp {crispr_result}/03.result/* {file_temp}" + "\n")
        w.write(f"rm {file_temp}/homo_region.result.xls" + "\n")
        w.write(f"head {crispr_result}/03.result/homo_region.result.xls >{file_temp}/homo_region.result.xls" + "\n")

def arrange_data_release(sh4):
    with open(sh4, "w") as w:
        w.write(f"cp -r {workflow_result}/published/data/01.fastq_qc {data_release}" + "\n")
        w.write(f"cp -r {workflow_result}/published/data/02.reference {data_release}" + "\n")
        w.write(f"cp -r {workflow_result}/published/data/03.mappingStat {data_release}" + "\n")
        w.write(f"cp -r {crispr_result}/03.result {data_release}/04.crispr_off_target" + "\n")
        w.write(f"cp -r {readme_dir}/04.crispr_off_target.readme.txt {data_release}/04.crispr_off_target/04.crispr_off_target.readme.txt" + "\n")

def generate_report(sh3):
    with open(sh3, "w") as w:
        w.write(f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report.rmd \
                --format html --outfile crispr_report.html")
        w.write("\n")
        w.write(f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report.rmd \
                --format pdf --outfile crispr_report.pdf")
        w.write("\n")
        w.write(f"cp {new_rmarkdown_scripts}/crispr_report.html {report_result}/{report_name}.html")
        w.write("\n")
        w.write(f"cp {new_rmarkdown_scripts}/crispr_report.pdf {report_result}/{report_name}.pdf")
        w.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="整理crispr_off_target的结果文件")
    parser.add_argument("-o", "--report_result", help="报告结果文件夹", required=True)
    parser.add_argument("-g", "--group_list", help="样本文件，便于获得数据", required=True)
    parser.add_argument("-w", "--workflow_result", help="变异检测的结果文件", required=True)
    parser.add_argument("-c", "--crispr_result", help="crispr脱靶率的结果文件", required=True)
    parser.add_argument("-r", "--rmarkdown_scripts", help="rmarkdown的路径",
                        default = os.path.join(os.path.split(os.path.realpath(__file__))[0], "rmarkdown_template"))
    parser.add_argument("-s", "--readme_dir", help="readme_dir路径",
                        default = os.path.join(os.path.split(os.path.realpath(__file__))[0], "readme_dir"))
    parser.add_argument("-n", "--report_name", help="报告名字", default="crispr_report")
    parser.add_argument("-a", "--only_report", help="是否只生成报告", action="store_true")
    argv=vars(parser.parse_args())
    '''
    全局配置
    '''
    rmarkdown_scripts = Path(argv["rmarkdown_scripts"])
    rmarkdown_scripts_bin = rmarkdown_scripts.joinpath("bin")
    group_list = argv["group_list"]
    workflow_result = argv["workflow_result"]
    crispr_result = argv["crispr_result"]
    readme_dir = argv["readme_dir"]
    report_name = argv["report_name"]
    """
    生成必要的文件
    """
    report_result = Path(argv["report_result"])
    report_result.mkdir(parents=True, exist_ok=True)
    new_rmarkdown_scripts = report_result.joinpath("rmarkdown_template")
    if not new_rmarkdown_scripts.exists():
        shutil.copytree(rmarkdown_scripts, new_rmarkdown_scripts)
    else:
        shutil.rmtree(new_rmarkdown_scripts)
        shutil.copytree(rmarkdown_scripts, new_rmarkdown_scripts)
    data_release = report_result.joinpath("data_release")
    data_release.mkdir(exist_ok=True)
    sh_temp = report_result.joinpath("sh_temp")
    sh_temp.mkdir(exist_ok=True)
    file_temp = report_result.joinpath("file")
    file_temp.mkdir(exist_ok=True)
    """
    wgs_v4变异检测结果整理
    """
    sh_1 = sh_temp.joinpath("1.sh")
    arrange_wgs_v4(sh_1)
    """
    脱靶率结果整理
    """
    sh_2 = sh_temp.joinpath("2.sh")
    arrange_crispr_off_target(sh_2)
    """
    结果文件整理
    """
    sh_4 = sh_temp.joinpath("4.sh")
    arrange_data_release(sh_4)
    """
    生成报告
    """
    sh_3 = sh_temp.joinpath("3.sh")
    generate_report(sh_3)
    """
    执行脚本
    """
    if argv["only_report"]:
        run_cmd(str(sh_3), f"只运行{sh_3}脚本，生产报告")
    else:
        sh_list = [sh_1, sh_2, sh_4, sh_3]
        for sh in sh_list:
            run_cmd(str(sh), f"运行{sh}脚本")
