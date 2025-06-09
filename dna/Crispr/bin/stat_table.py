# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/08/25
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import csv
from collections import namedtuple


def stat_table(infile, outfile):
    #infile="/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_development/crispr_off_target/20230818_test/demo/02.filter_vcf/targe_info.xls"
    site_list = []
    site_snv_list = []
    site_indel_list = []
    gene_list = []
    with open(infile, "r") as f:
        f_csv = csv.reader(f, delimiter="\t")
        header = next(f_csv)
        header = [i.replace(".", "_") for i in header]
        table_info_nt = namedtuple("table_info_nt", header)
        for row in f_csv:
            each_row = table_info_nt(*row)
            site = each_row.CHROM + each_row.POS
            if each_row.Vtype == "SNV":
                site_snv_list.append(site)
            else:
                site_indel_list.append(site)
            gene = each_row.GeneName
            site_list.append(site)
            gene_list.append(gene)
    site_num = len(set(site_list))
    snv_site_num = len(set(site_snv_list))
    indel_site_num = len(set(site_indel_list))
    gene_num = len(set(gene_list))
    with open(outfile, "w") as w:
        header_list = ["Total potential off-targets", "Potential SNPs", "Potential InDels", "Genes"]
        line_list = [str(site_num), str(snv_site_num), str(indel_site_num), str(gene_num)]
        w.write("\t".join(header_list) + "\n")
        w.write("\t".join(line_list) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="table统计")
    parser.add_argument("-i", "--infile", help="table表", required=True)
    parser.add_argument("-o", "--outfile", help="输出文件", required=True)
    args = parser.parse_args()
    stat_table(args.infile, args.outfile)
