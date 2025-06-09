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


def arrange_vcf2table(sh1, src_dir, dest_dir):
    with open(sh1, "w") as w:
        w.write(
            """awk '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5}' %s/pop.index | cut -f 3-5 > %s/tmp2.txt"""
            % (src_dir, dest_dir)
        )
        w.write("\n")
        w.write(f"cat {src_dir}/pop.table |cut -f 1-5 > {dest_dir}/tmp1.txt")
        w.write("\n")
        w.write(
            """awk -F "\\t" '{for (i=6;i<=NF;i++)printf("%%s\\t", $i);print ""}' %s/pop.table > %s/tmp3.txt"""
            % (src_dir, dest_dir)
        )
        w.write("\n")
        w.write(
            f"paste {dest_dir}/tmp1.txt {dest_dir}/tmp2.txt {dest_dir}/tmp3.txt > {dest_dir}/pop.table"
        )
        w.write("\n")
        w.write(f"rm {dest_dir}/tmp1.txt {dest_dir}/tmp2.txt {dest_dir}/tmp3.txt")
        w.write("\n")
        w.write(f"cp {src_dir}/snp_indel_gene.stat {dest_dir}/snp_indel_gene.stat.xls")
        w.write("\n")
        w.write(
            f"""sed -i '1s/.*/CHROM\\tSNP Number\\tEffective SNP\\tInDel Number\\tEffective InDel\\tGene Number\\tEffective Gene/' {dest_dir}/snp_indel_gene.stat.xls"""
        )
        w.write("\n")
        w.write(f"cp {src_dir}/pop.final.anno {dest_dir}/pop.final.anno")
        w.write("\n")
        w.write(
            f"Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {dest_dir}/pop.final.anno --outfile {dest_dir}/pop.final.anno.xls"
        )
        w.write("\n")
        w.write(f"rm -rf {dest_dir}/pop.final.anno")
        w.write("\n")
        w.write(f"rm -rf {dest_dir}/pop.final.anno")
        w.write("\n")


def arrange_wgs2table(sh17, src_dir, dest_dir):
    with open(sh17, "w") as w:
        w.write(f"mkdir -p {dest_dir}/region")
        w.write("\n")
        w.write(f"cp {src_dir}/*.xls {dest_dir}/region/")
        w.write("\n")


def arrange_index_slid(sh2, src_dir, dest_dir):
    with open(sh2, "w") as w:
        w.write(
            """awk '{print $1"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}' %s/pop.bootstrap.result > %s/pop.index.bootstrap.result.xls"""
            % (src_dir, dest_dir)
        )
        w.write("\n")
        w.write(
            f"""for i in $(cut -f1 {chr_list});do
                cp {src_dir}/pop.index.$i.index.png {dest_dir}/index.$i.png
                cp {src_dir}/pop.index.$i.index.pdf {dest_dir}/index.$i.pdf
                done
                """
        )
        w.write("\n")
        w.write(f"cp {src_dir}/pop.index.index.png {dest_dir}/index.png")
        w.write("\n")
        w.write(f"cp {src_dir}/pop.index.index.pdf {dest_dir}/index.pdf")
        w.write("\n")
        w.write(
            f"""if [ -d {enrich_result}/index ];then
                    Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/index/index.all.table --outfile {dest_dir}/index.all.table.xls
                    cp {enrich_result}/index/index.region_stat.xls {dest_dir}/index.region_stat.xls
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""sed "s/ /\\t/g" {src_dir}/pop.bootstrap.result > {dest_dir}/index.detail.result.xls"""
        )
        w.write("\n")


def arrange_loess(sh3, src_dir, dest_dir):
    with open(sh3, "w") as w:
        w.write(
            """awk '{print $1"\\t"$2"\\t"$8"\\t"$9"\\t"$10}' %s/pop.bootstrap.result > %s/pop.loess.bootstrap.result.xls"""
            % (src_dir, dest_dir)
        )
        w.write("\n")
        w.write(
            f"""for i in $(cut -f1 {chr_list});do
                cp {src_dir}/loess.$i.index.png {dest_dir}/loess.$i.png
                cp {src_dir}/loess.$i.index.pdf {dest_dir}/loess.$i.pdf
                done
                """
        )
        w.write("\n")
        w.write(f"cp {src_dir}/loess.index.png {dest_dir}/loess.png")
        w.write("\n")
        w.write(f"cp {src_dir}/loess.index.pdf {dest_dir}/loess.pdf")
        w.write("\n")
        w.write(
            f"""if [ -d {enrich_result}/loess ];then
                    Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/loess/loess.all.table --outfile {dest_dir}/loess.all.table.xls
                    cp {enrich_result}/loess/loess.region_stat.xls {dest_dir}/loess.region_stat.xls
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""sed "s/ /\\t/g" {src_dir}/pop.bootstrap.result > {dest_dir}/loess.detail.result.xls"""
        )
        w.write("\n")


def arrange_ed(sh4, src_dir, dest_dir, mutmap):
    with open(sh4, "w") as w:
        if mutmap:
            w.write("""echo Mutmap BSA. No ED.""")
        else:
            w.write(
                """awk '{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$8}' %s/pop.sliding.result > %s/pop.ED.sliding.result.xls"""
                % (src_dir, dest_dir)
            )
            w.write("\n")
            w.write(
                f"""for i in $(cut -f1 {chr_list});do
                    cp {src_dir}/pop.ED.$i.index.png {dest_dir}/ED.$i.png
                    cp {src_dir}/pop.ED.$i.index.pdf {dest_dir}/ED.$i.pdf
                    done
                    """
            )
            w.write("\n")
            w.write(f"cp {src_dir}/pop.ED.index.png {dest_dir}/ED.png")
            w.write("\n")
            w.write(f"cp {src_dir}/pop.ED.index.pdf {dest_dir}/ED.pdf")
            w.write("\n")
            w.write(
                f"""if [ -d {enrich_result}/ED ];then
                        Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/ED/ED.all.table --outfile {dest_dir}/ED.all.table.xls
                        cp {enrich_result}/ED/ED.region_stat.xls {dest_dir}/ED.region_stat.xls
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""sed "s/ /\\t/g" {src_dir}/pop.sliding.detail | cut -f1,2,13,14,16,17,18,19,20,21,22,23,24,25 > {dest_dir}/ED.detail.result.xls"""
            )
            w.write("\n")


def arrange_gprime(sh5, src_dir, dest_dir, mutmap):
    with open(sh5, "w") as w:
        if mutmap:
            w.write("""echo Mutmap BSA. No Gprime.""")
        else:
            w.write(
                """awk '{print $1"\\t"$2"\\t"$12"\\t"$13"\\t"$15"\\t"$16"\\t"$17}' %s/pop.Gprime > %s/pop.Gprime.result.xls"""
                % (src_dir, dest_dir)
            )
            w.write("\n")
            w.write(
                f"""for i in $(cut -f1 {chr_list});do
                    cp {src_dir}/Gprime.$i.index.png {dest_dir}/Gprime.$i.png
                    cp {src_dir}/Gprime.$i.index.pdf {dest_dir}/Gprime.$i.pdf
                    done
                    """
            )
            w.write("\n")
            w.write(f"cp {src_dir}/Gprime.index.png {dest_dir}/Gprime.png")
            w.write("\n")
            w.write(f"cp {src_dir}/Gprime.index.pdf {dest_dir}/Gprime.pdf")
            w.write("\n")
            w.write(
                f"""if [ -d {enrich_result}/Gprime ];then
                        Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/Gprime/Gprime.all.table --outfile {dest_dir}/Gprime.all.table.xls
                        cp {enrich_result}/Gprime/Gprime.region_stat.xls {dest_dir}/Gprime.region_stat.xls
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""sed "s/ /\\t/g" {src_dir}/pop.Gprime > {dest_dir}/Gprime.detail.result.xls"""
            )
            w.write("\n")


def arrange_deepbsa_dl(sh6, src_dir, dest_dir, include_deepbsa):
    with open(sh6, "w") as w:
        if include_deepbsa:
            w.write(f"cp {src_dir}/DL_values.txt {dest_dir}/pop.DeepBSA_DL.result.xls")
            w.write("\n")
            w.write(
                f"""for i in $(cut -f1 {chr_list});do
                    cp {src_dir}/pop.deepBSA_DL.$i.index.png {dest_dir}/DeepBSA_DL.$i.png
                    cp {src_dir}/pop.deepBSA_DL.$i.index.pdf {dest_dir}/DeepBSA_DL.$i.pdf
                    done
                    """
            )
            w.write("\n")
            w.write(f"cp {src_dir}/pop.deepBSA_DL.index.png {dest_dir}/DeepBSA_DL.png")
            w.write("\n")
            w.write(f"cp {src_dir}/pop.deepBSA_DL.index.pdf {dest_dir}/DeepBSA_DL.pdf")
            w.write("\n")
            w.write(
                f"""if [ -d {enrich_result}/DeepBSA_DL ];then
                        Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/DeepBSA_DL/DeepBSA_DL.all.table --outfile {dest_dir}/DeepBSA_DL.all.table.xls
                        cp {enrich_result}/DeepBSA_DL/DeepBSA_DL.region_stat.xls {dest_dir}/DeepBSA_DL.region_stat.xls
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""sed "s/ /\\t/g" {src_dir}/DL_values.txt > {dest_dir}/DeepBSA_DL.detail.result.xls"""
            )
            w.write("\n")
        else:
            w.write("echo No DeepBSA.")


def arrange_deepbsa_k(sh7, src_dir, dest_dir, include_deepbsa):
    with open(sh7, "w") as w:
        if include_deepbsa:
            w.write(f"cp {src_dir}/K_values.txt {dest_dir}/pop.DeepBSA_K.result.xls")
            w.write("\n")
            w.write(
                f"""for i in $(cut -f1 {chr_list});do
                    cp {src_dir}/pop.deepBSA_K.$i.index.png {dest_dir}/DeepBSA_K.$i.png
                    cp {src_dir}/pop.deepBSA_K.$i.index.pdf {dest_dir}/DeepBSA_K.$i.pdf
                    done
                    """
            )
            w.write("\n")
            w.write(f"cp {src_dir}/pop.deepBSA_K.index.png {dest_dir}/DeepBSA_K.png")
            w.write("\n")
            w.write(f"cp {src_dir}/pop.deepBSA_K.index.pdf  {dest_dir}/DeepBSA_K.pdf")
            w.write("\n")
            w.write(
                f"""if [ -d {enrich_result}/DeepBSA_DL ];then
                        Rscript {rmarkdown_scripts_bin}/process_final_anno.r --infile {enrich_result}/DeepBSA_K/DeepBSA_K.all.table --outfile {dest_dir}/DeepBSA_K.all.table.xls
                        cp {enrich_result}/DeepBSA_K/DeepBSA_K.region_stat.xls {dest_dir}/DeepBSA_K.region_stat.xls
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""sed "s/ /\\t/g" {src_dir}/K_values.txt > {dest_dir}/DeepBSA_K.detail.result.xls"""
            )
            w.write("\n")
        else:
            w.write("echo No DeepBSA.")


def arrange_enrich(sh8, src_dir, dest_dir, include_deepbsa, mutmap):
    with open(sh8, "w") as w:
        if mutmap:
            w.write("for i in index loess;do")
        elif include_deepbsa:
            w.write("for i in ED Gprime index loess DeepBSA_DL DeepBSA_K;do")
        else:
            w.write("for i in ED Gprime index loess;do")

        w.write("\n")
        w.write(
            f"""mkdir -p {dest_dir}/$i
                    cp -r {src_dir}/$i/GO_result {dest_dir}/$i/
                    cp -r {src_dir}/$i/KEGG_result {dest_dir}/$i/
                    cp {src_dir}/$i.degfile {dest_dir}/$i/$i.degfile.xls
                    cp {src_dir}/$i.transcript.degfile {dest_dir}/$i/$i.transcript.degfile.xls
                    cp {src_dir}/$i.genes_abstract.list {dest_dir}/$i/$i.genes_abstract.list.xls
                    cp {src_dir}/$i.extract.raw.vcf.gz {dest_dir}/$i/$i.extract.raw.vcf.gz
                    done
                """
        )
        w.write("\n")


def arrange_file(sh9, src_dir, dest_dir, deepbsa, mutmap):
    with open(sh9, "w") as w:
        w.write(f"cp {group_list} {dest_dir}/group.list")
        w.write("\n")
        w.write(f"cp {chr_list} {dest_dir}/chr.list")
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/02.index/index.region_stat.xls ];then
                   sed '1d' {src_dir}/02.index/index.region_stat.xls | sed 's/^/index-slid\\t/g' > {dest_dir}/index.region_stat
                   fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/03.loess/loess.region_stat.xls ];then
                    sed '1d' {src_dir}/03.loess/loess.region_stat.xls | sed 's/^/index-loess\\t/g' > {dest_dir}/loess.region_stat
                    fi
                """
        )
        w.write("\n")
        if not mutmap:
            w.write(
                f"""if [ -e {src_dir}/04.ED/ED.region_stat.xls ];then
                        sed '1d' {src_dir}/04.ED/ED.region_stat.xls | sed 's/^/Euclidean\\t/g' > {dest_dir}/ED.region_stat
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""if [ -e {src_dir}/05.Gprime/Gprime.region_stat.xls ];then
                        sed '1d' {src_dir}/05.Gprime/Gprime.region_stat.xls | sed 's/^/Gprime\\t/g' > {dest_dir}/Gprime.region_stat
                        fi
                    """
            )
            w.write("\n")
        if deepbsa:
            w.write(
                f"""if [ -e {src_dir}/06.DeepBSA_DL/DeepBSA_DL.region_stat.xls ];then
                        sed '1d' {src_dir}/06.DeepBSA_DL/DeepBSA_DL.region_stat.xls | sed 's/^/DeepBSA_DL\\t/g' > {dest_dir}/DeepBSA_DL.region_stat
                        fi
                    """
            )
            w.write("\n")
            w.write(
                f"""if [ -e {src_dir}/07.DeepBSA_K/DeepBSA_K.region_stat.xls ];then
                        sed '1d' {src_dir}/07.DeepBSA_K/DeepBSA_K.region_stat.xls | sed 's/^/DeepBSA_K\\t/g' > {dest_dir}/DeepBSA_K.region_stat
                        fi
                    """
            )
            w.write("\n")
        w.write(
            f"""cat {dest_dir}/*.region_stat | sed '1i\Method\\tRegion\\tSNP_Number\\tEffective_SNP\\tInDel_Number\\tEffective_InDel\\tGene_Number\\tNRID\\tUniID\\tKoID\\tGOTERM\\tEGGNOG\\tPfamID'> {dest_dir}/all.region_stat.xls"""
        )
        w.write("\n")
        w.write(
            f"cp {src_dir}/01.vcf2table/snp_indel_gene.stat.xls {dest_dir}/snp_indel_gene.stat.xls"
        )
        w.write("\n")
        w.write(f"cp {src_dir}/02.index/index.png {dest_dir}")
        w.write("\n")
        w.write(f"cp {src_dir}/03.loess/loess.png {dest_dir}")
        w.write("\n")
        if not mutmap:
            w.write(f"cp {src_dir}/04.ED/ED.png {dest_dir}")
            w.write("\n")
            w.write(f"cp {src_dir}/05.Gprime/Gprime.png {dest_dir}")
            w.write("\n")
        if deepbsa:
            w.write(f"cp {src_dir}/06.DeepBSA_DL/DeepBSA_DL.png {dest_dir}")
            w.write("\n")
            w.write(f"cp {src_dir}/07.DeepBSA_K/DeepBSA_K.png {dest_dir}")
            w.write("\n")


def arrange_file_enrich(sh10, src_dir, dest_dir):
    with open(sh10, "w") as w:
        w.write(
            f"""if [ -e {src_dir}/index/GO_result/index_GOenrichment.png ];then
                    cp {src_dir}/index/GO_result/index_GOenrichment.png {dest_dir}/index_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/loess/GO_result/loess_GOenrichment.png ];then
                    cp {src_dir}/loess/GO_result/loess_GOenrichment.png {dest_dir}/loess_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/ED/GO_result/ED_GOenrichment.png ];then
                    cp {src_dir}/ED/GO_result/ED_GOenrichment.png {dest_dir}/ED_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/Gprime/GO_result/Gprime_GOenrichment.png ];then
                    cp {src_dir}/Gprime/GO_result/Gprime_GOenrichment.png {dest_dir}/Gprime_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/DeepBSA_DL/GO_result/DeepBSA_DL_GOenrichment.png ];then
                    cp {src_dir}/DeepBSA_DL/GO_result/DeepBSA_DL_GOenrichment.png {dest_dir}/DeepBSA_DL_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/DeepBSA_K/GO_result/DeepBSA_K_GOenrichment.png ];then
                    cp {src_dir}/DeepBSA_K/GO_result/DeepBSA_K_GOenrichment.png {dest_dir}/DeepBSA_K_go_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/index/KEGG_result/index_KEGGenrichment.png ];then
                    cp {src_dir}/index/KEGG_result/index_KEGGenrichment.png {dest_dir}/index_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/loess/KEGG_result/loess_KEGGenrichment.png ];then
                    cp {src_dir}/loess/KEGG_result/loess_KEGGenrichment.png {dest_dir}/loess_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/ED/KEGG_result/ED_KEGGenrichment.png ];then
                    cp {src_dir}/ED/KEGG_result/ED_KEGGenrichment.png {dest_dir}/ED_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/Gprime/KEGG_result/Gprime_KEGGenrichment.png ];then
                    cp {src_dir}/Gprime/KEGG_result/Gprime_KEGGenrichment.png {dest_dir}/Gprime_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/DeepBSA_DL/KEGG_result/DeepBSA_DL_KEGGenrichment.png ];then
                    cp {src_dir}/DeepBSA_DL/KEGG_result/DeepBSA_DL_KEGGenrichment.png {dest_dir}/DeepBSA_DL_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(
            f"""if [ -e {src_dir}/DeepBSA_K/KEGG_result/DeepBSA_K_KEGGenrichment.png ];then
                    cp {src_dir}/DeepBSA_K/KEGG_result/DeepBSA_K_KEGGenrichment.png {dest_dir}/DeepBSA_K_kegg_enrich.png
                    fi
                """
        )
        w.write("\n")
        w.write(f"cp {src_dir}/*/*.degfile.xls {dest_dir}")


def arrange_wgs_v4(sh11, src_dir, dest_dir):
    with open(sh11, "w") as w:
        ## qc结果整理
        w.write(
            f"""for i in `cut -f 1 {group_list}`
                    do
                    cp {src_dir}/published/data/01.fastq_qc/fig/$i.raw.base.png {dest_dir}
                    cp {src_dir}/published/data/01.fastq_qc/fig/$i.clean.base.png {dest_dir}
                    cp {src_dir}/published/data/01.fastq_qc/fig/$i.raw.qual.png {dest_dir}
                    cp {src_dir}/published/data/01.fastq_qc/fig/$i.clean.qual.png {dest_dir}
                    done
                """
        )
        w.write("\n")
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_qc_result.py --qc_stat {src_dir}/published/data/01.fastq_qc/qc.stat.txt \
                --group_file {group_list} --raw_data {dest_dir}/rawdata.xls --clean_data {dest_dir}/cleandata.xls"
        )
        w.write("\n")
        ## 比对结果统计
        w.write(
            f"""for i in `cut -f 1 {group_list}`
                    do
                    cp {src_dir}/published/data/03.mappingStat/$i.genome.coverage.png {dest_dir}
                    cp {src_dir}/published/data/03.mappingStat/$i.insert.png {dest_dir}
                    cp {src_dir}/published/data/03.mappingStat/$i.depth.png {dest_dir}
                    done
                """
        )
        w.write("\n")
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_mapping_result.py --all_summary_stats {src_dir}/published/data/03.mappingStat/all.summary.stats.xls \
                --group_file {group_list} --bsa_align_stat {dest_dir}/align_stat.xls"
        )
        w.write("\n")
        ## 变异检测结果整理
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_snp_result.py --group_file {group_list} \
                --snp_stat {src_dir}/published/data/04.snpIndel/snp/snp.stat.xls --snp_stat_bsa {dest_dir}/snp_stat.xls"
        )
        w.write("\n")
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_indel_result.py --group_file {group_list} \
                --indel_stat {src_dir}/published/data/04.snpIndel/indel/indel.stat.xls --indel_stat_bsa {dest_dir}/indel_stat.xls"
        )
        w.write("\n")
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_snp_anno_result.py --group_file {group_list} \
                --snp_anno_stat {src_dir}/published/data/04.snpIndel/snp/snp_anno.stat.xls \
                --snp_anno_stat_bsa {dest_dir}/snp_anno.xls"
        )
        w.write("\n")
        w.write(
            """cut -f 1,10,12,17-19,21 %s/snp_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
                | sed 's/[ \\t]*$//g' | sed 's/ /\\t/g' > %s/snp_anno_T.xls"""
            % (dest_dir, dest_dir)
        )
        w.write("\n")
        w.write(f"cut -f 1,5-8 {dest_dir}/snp_anno.xls > {dest_dir}/snp_impact.xls")
        w.write("\n")
        w.write(
            f"python3 {rmarkdown_scripts_bin}/process_indel_anno_result.py --group_file {group_list} \
                --indel_anno_stat {src_dir}/published/data/04.snpIndel/indel/indel_anno.stat.xls \
                --indel_anno_stat_bsa {dest_dir}/indel_anno.xls"
        )
        w.write("\n")
        w.write(
            """cut -f 1,16,17,19,27,29,30 %s/indel_anno.xls | awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' \
                | sed 's/[ \\t]*$//g' | sed 's/ /\\t/g' > %s/indel_anno_T.xls"""
            % (dest_dir, dest_dir)
        )
        w.write("\n")
        w.write(f"cut -f 1,6-9 {dest_dir}/indel_anno.xls > {dest_dir}/indel_impact.xls")
        w.write("\n")
        w.write(f"cp {src_dir}/tmp/02.reference/project.info {dest_dir}")
        w.write("\n")
        w.write(f"cp {src_dir}/tmp/02.reference/info.log {dest_dir}")
        w.write("\n")


def arrange_wgs_v4_other(sh12, src_dir, dest_dir):
    with open(sh12, "w") as w:
        w.write(f"cp {src_dir}/snp_anno.xls {dest_dir}")
        w.write("\n")
        w.write(f"cp {src_dir}/indel_anno.xls {dest_dir}")
        w.write("\n")
        w.write(
            f"cp {workflow_result}/published/data/02.reference/ref.genome.summary.xls {dest_dir}"
        )


def generate_report(sh13, dest_dir, report_name, rna: bool, mutmap, deepbsa):
    with open(sh13, "w") as w:
        if rna:
            w.write(
                f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report_rna.rmd \
                    --format html --outfile {new_rmarkdown_scripts}/bsa_report.html --population {population} --primer_design {primer_design} {'--mutmap yes' if mutmap else '--mutmap no'} {'--deepbsa yes' if deepbsa else '--deepbsa no'}"
            )
            w.write("\n")
            w.write(
                f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report_rna.rmd \
                    --format pdf --outfile {new_rmarkdown_scripts}/bsa_report.pdf --population {population} --primer_design {primer_design} {'--mutmap yes' if mutmap else '--mutmap no'} {'--deepbsa yes' if deepbsa else '--deepbsa no'}"
            )
            w.write("\n")
        else:
            w.write(
                f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report.rmd \
                    --format html --outfile {new_rmarkdown_scripts}/bsa_report.html --population {population} --primer_design {primer_design} {'--mutmap yes' if mutmap else '--mutmap no'} {'--deepbsa yes' if deepbsa else '--deepbsa no'}"
            )
            w.write("\n")
            w.write(
                f"Rscript {new_rmarkdown_scripts}/rmarkdown.r --rmd  {new_rmarkdown_scripts}/report.rmd \
                    --format pdf --outfile {new_rmarkdown_scripts}/bsa_report.pdf --population {population} --primer_design {primer_design} {'--mutmap yes' if mutmap else '--mutmap no'} {'--deepbsa yes' if deepbsa else '--deepbsa no'}"
            )
            w.write("\n")
        w.write(
            f"cp {new_rmarkdown_scripts}/bsa_report.html {dest_dir}/{report_name}.html"
        )
        w.write("\n")
        w.write(
            f"cp {new_rmarkdown_scripts}/bsa_report.pdf {dest_dir}/{report_name}.pdf"
        )
        w.write("\n")


def generate_readme(sh14, dest_dir, include_deepbsa):
    with open(sh14, "w") as w:
        w.write(f"cp -r {new_rmarkdown_scripts}/readme_dir_qtl {dest_dir}")
        w.write("\n")
        if include_deepbsa:
            w.write(
                f"bash {new_rmarkdown_scripts}/generate_readme_use_deepbsa.bash {dest_dir} {chr_list}"
            )
            w.write("\n")
            w.write(
                f"python3 {dest_dir}/readme_dir_qtl/ReadMe_result.py --file {dest_dir}/BSA.info.result --typeinfo BSA_new --outname {data_release}/{report_name} --readme_dir {new_rmarkdown_scripts}/readme_dir_qtl --data_release {data_release}"
            )
        else:
            w.write(
                f"bash {new_rmarkdown_scripts}/generate_readme_no_use_deepbsa.bash {dest_dir} {chr_list}"
            )
            w.write("\n")
            w.write(
                f"python3 {dest_dir}/readme_dir_qtl/ReadMe_result.py --file {dest_dir}/BSA.info.result --typeinfo BSA_new_no_use_deepbsa --outname {data_release}/结果目录索引及说明 --readme_dir {new_rmarkdown_scripts}/readme_dir_qtl  --data_release {data_release}"
            )
            w.write("\n")


def arrange_primer_design(sh15, src_dir, include_deepbsa, primer_design):
    if primer_design:
        with open(sh15, "w") as w:
            if include_deepbsa:
                w.write(
                    f"mkdir -p {data_release}/09.primer_design"
                )
                w.write("\n")
                w.write(
                    f"cp -r {primer_design_result}/07.result/* {data_release}/09.primer_design"
                )
                w.write("\n")
                w.write(
                    f"python3 {rmarkdown_scripts_bin}/stat_primer_result.py --indir {data_release}/09.primer_design --outfile {file_temp}/primer_design_stat.xls"
                )
            else:
                w.write(
                    f"mkdir -p {data_release}/07.primer_design"
                )
                w.write("\n")
                w.write(
                    f"cp -r {primer_design_result}/07.result/* {data_release}/07.primer_design"
                )
                w.write("\n")
                w.write(
                    f"python3 {rmarkdown_scripts_bin}/stat_primer_result.py --indir {data_release}/07.primer_design --outfile {file_temp}/primer_design_stat.xls"
                )
    else:
        return


def arrange_params_info(sh16):
    with open(sh16, "w") as w:
        w.write(
            f"cat {bsa_result}/*/*.params.log  > {file_temp}/params.xls"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="整理BSA的结果文件")
    parser.add_argument("-i", "--bsa_result_dir", help="输入文件", required=True)
    parser.add_argument("-o", "--report_result", help="报告结果文件夹", required=True)
    parser.add_argument(
        "-e", "--include_deepBSA", help="是否不带deepBSA方法", action="store_true"
    )
    parser.add_argument("-g", "--group_list", help="bsa的分组文件", required=True)
    parser.add_argument("-c", "--chr_list", help="要作图的chr_list文件", required=True)
    parser.add_argument("-w", "--workflow_result", help="变异检测的结果文件", required=True)
    parser.add_argument(
        "-r",
        "--rmarkdown_scripts",
        help="rmarkdown的路径",
        default=os.path.join(
            os.path.split(os.path.realpath(__file__))[0], "rmarkdown_template_pdf"
        ),
    )
    parser.add_argument("-n", "--report_name", help="报告名字", default="BSA分析报告")
    parser.add_argument("-p", "--population", help="群体类型", default="F2")
    parser.add_argument("-a", "--only_report", help="是否只生成报告", action="store_true")
    parser.add_argument("-b", "--only_readme", help="是否只生产readme", action="store_true")
    parser.add_argument(
        "-d", "--primer_design", help="是否包含引物设计的报告", action="store_true"
    )
    parser.add_argument("-j", "--primer_design_result", help="引物设计结果")
    parser.add_argument("-k", "--rna", help="BSR报告", action="store_true", default=False)
    parser.add_argument(
        "-m", "--mutmap", help="mutmap报告", action="store_true", default=False
    )
    argv = vars(parser.parse_args())
    """
    全局配置
    """
    population = argv["population"]
    if argv["primer_design"]:
        primer_design = "yes"
        primer_design_result = Path(argv["primer_design_result"]).absolute()
    else:
        primer_design = "no"
    rmarkdown_scripts = Path(argv["rmarkdown_scripts"]).absolute()
    rmarkdown_scripts_bin = rmarkdown_scripts.joinpath("bin")
    chr_list = str(Path(argv["chr_list"]).absolute())
    group_list = str(Path(argv["group_list"]).absolute())
    workflow_result = str(Path(argv["workflow_result"]).absolute())
    report_name = argv["report_name"]
    if argv["rna"] and argv["report_name"] == "BSA分析报告":
        report_name = "BSR分析报告"
    """
    生成必要的文件
    """
    bsa_result = Path(argv["bsa_result_dir"]).absolute()
    report_result = Path(argv["report_result"]).absolute()
    if not bsa_result.exists():
        raise Exception("BSA结果文件夹不存在")
    report_result.mkdir(parents=True, exist_ok=True)
    new_rmarkdown_scripts = report_result.joinpath("rmarkdown_template_pdf")
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
    if argv["include_deepBSA"]:
        release_dir_list = [
            "01.vcf2table",
            "02.index",
            "03.loess",
            "04.ED",
            "05.Gprime",
            "06.DeepBSA_DL",
            "07.DeepBSA_K",
            "08.enrich",
        ]
        enrich_result = bsa_result.joinpath("06.enrich")
        src_vcf2table = bsa_result.joinpath("01.vcf2table")
        dest_vcf2table = data_release.joinpath("01.vcf2table")
        src_index = bsa_result.joinpath("02.index-slid")
        dest_index = data_release.joinpath("02.index")
        src_loess = bsa_result.joinpath("05.loess")
        dest_loess = data_release.joinpath("03.loess")
        src_ed = bsa_result.joinpath("04.ED-slid")
        dest_ed = data_release.joinpath("04.ED")
        src_gprime = bsa_result.joinpath("03.Gprime")
        dest_gprime = data_release.joinpath("05.Gprime")
        src_deepbsa = bsa_result.joinpath("08.DeepBSA")
        dest_deepbsa_dl = data_release.joinpath("06.DeepBSA_DL")
        dest_deepbsa_k = data_release.joinpath("07.DeepBSA_K")
        dest_enrich = data_release.joinpath("08.enrich")
        src_primer_design = bsa_result.joinpath("09.primer_design")
        # 新增region区域增加anno注释，gene注释
        wgs_vcf2table = bsa_result.joinpath("07.extract_vcf")
    else:
        release_dir_list = [
            "01.vcf2table",
            "02.index",
            "03.loess",
            "04.ED",
            "05.Gprime",
            "06.enrich",
        ]
        enrich_result = bsa_result.joinpath("06.enrich")
        src_vcf2table = bsa_result.joinpath("01.vcf2table")
        dest_vcf2table = data_release.joinpath("01.vcf2table")
        src_index = bsa_result.joinpath("02.index-slid")
        dest_index = data_release.joinpath("02.index")
        src_loess = bsa_result.joinpath("05.loess")
        dest_loess = data_release.joinpath("03.loess")
        src_ed = bsa_result.joinpath("04.ED-slid")
        dest_ed = data_release.joinpath("04.ED")
        src_gprime = bsa_result.joinpath("03.Gprime")
        dest_gprime = data_release.joinpath("05.Gprime")
        dest_enrich = data_release.joinpath("06.enrich")
        src_deepbsa = ""
        dest_deepbsa_dl = ""
        dest_deepbsa_k = ""
        src_primer_design = bsa_result.joinpath("09.primer_design")
        # 新增region区域增加anno注释，gene注释
        wgs_vcf2table = bsa_result.joinpath("07.extract_vcf")
    for dir in release_dir_list:
        data_release.joinpath(dir).mkdir(exist_ok=True)
    """
    整理01.vcf2table文件夹
    """
    sh_1 = sh_temp.joinpath("1.sh")
    arrange_vcf2table(sh_1, src_vcf2table, dest_vcf2table)
    """
    整理01.vcf2table文件夹,增加region区域
    """
    sh_17 = sh_temp.joinpath("17.sh")
    arrange_wgs2table(sh_17, wgs_vcf2table, dest_enrich)
    """
    整理02.index文件夹
    """
    sh_2 = sh_temp.joinpath("2.sh")
    arrange_index_slid(sh_2, src_index, dest_index)
    """
    整理03.loess文件夹
    """
    sh_3 = sh_temp.joinpath("3.sh")
    arrange_loess(sh_3, src_loess, dest_loess)
    """
    整理04.ED文件夹
    """
    sh_4 = sh_temp.joinpath("4.sh")
    arrange_ed(sh_4, src_ed, dest_ed, argv["mutmap"])
    """
    整理05.Gprime结果整理
    """
    sh_5 = sh_temp.joinpath("5.sh")
    arrange_gprime(sh_5, src_gprime, dest_gprime, argv["mutmap"])
    """
    整理06.DeepBSA_DL文件夹
    """
    sh_6 = sh_temp.joinpath("6.sh")
    arrange_deepbsa_dl(sh_6, src_deepbsa, dest_deepbsa_dl, argv["include_deepBSA"])
    """
    整理07.DeepBSA_K文件夹
    """
    sh_7 = sh_temp.joinpath("7.sh")
    arrange_deepbsa_k(sh_7, src_deepbsa, dest_deepbsa_k, argv["include_deepBSA"])
    """
    整理08.enrich文件夹
    """
    sh_8 = sh_temp.joinpath("8.sh")
    arrange_enrich(sh_8, enrich_result, dest_enrich, argv["include_deepBSA"], argv["mutmap"])
    """
    整理文件至file文件夹,方便rmd调用
    """
    sh_9 = sh_temp.joinpath("9.sh")
    arrange_file(sh_9, data_release, file_temp, argv["include_deepBSA"], argv["mutmap"])
    """
    整理富集文件图到file文件夹,方便rmd调用
    """
    sh_10 = sh_temp.joinpath("10.sh")
    arrange_file_enrich(sh_10, dest_enrich, file_temp)
    """
    wgs_v4变异检测结果整理
    """
    sh_11 = sh_temp.joinpath("11.sh")
    arrange_wgs_v4(sh_11, workflow_result, file_temp)
    """
    整理project.info和info.log
    """
    sh_12 = sh_temp.joinpath("12.sh")
    arrange_wgs_v4_other(sh_12, file_temp, dest_vcf2table)
    """
    生成报告
    """
    sh_13 = sh_temp.joinpath("13.sh")
    generate_report(
        sh_13,
        report_result,
        report_name,
        argv["rna"],
        argv["mutmap"],
        argv["include_deepBSA"]
    )
    """
    生成readme脚本
    """
    sh_14 = sh_temp.joinpath("14.sh")
    generate_readme(sh_14, report_result, argv["include_deepBSA"])
    """
    整理引物设计文件夹
    """
    sh_15 = sh_temp.joinpath("15.sh")
    arrange_primer_design(
        sh_15, src_primer_design, argv["include_deepBSA"], argv["primer_design"]
    )
    sh_16 = sh_temp.joinpath("16.sh")
    arrange_params_info(sh_16)
    """
    执行脚本
    """
    if argv["only_report"]:
        run_cmd(str(sh_13), f"只运行{sh_13}脚本，生产报告")
    if argv["only_readme"]:
        run_cmd(str(sh_14), f"只运行{sh_14}脚本，生产报告")
    else:
        sh_list = [
            sh_1,
            sh_2,
            sh_3,
            sh_4,
            sh_5,
            sh_6,
            sh_7,
            sh_8,
            sh_9,
            sh_10,
            sh_15,
            sh_16,
            sh_11,
            sh_12,
            sh_13,
            sh_14,
            sh_17,
        ]
        for sh in sh_list:
            run_cmd(str(sh), f"运行{sh}脚本")
