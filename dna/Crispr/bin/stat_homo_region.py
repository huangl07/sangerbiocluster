# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/08/31
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""


import argparse
import csv
from collections import namedtuple, defaultdict


def stat_homo_region(infile, outfile):
    #infile="/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_development/crispr_off_target/20230831_demo_test/demo/03.result/homo_region.result.xls"
    site_dict = defaultdict(list)
    with open(infile, "r") as f:
        f_csv = csv.reader(f, delimiter="\t")
        header = next(f_csv)
        header = [i.replace(" ", "_") for i in header]
        homo_region_nt = namedtuple("homo_region_nt", header)
        for row in f_csv:
            each_row = homo_region_nt(*row)
            mismatch = each_row.Mismatches
            bulge_size = each_row.Bulge_Size
            chrom = each_row.Chromosome
            position = each_row.Position
            key = (bulge_size, mismatch)
            site_dict[key].append(chrom + "_" + position)
    sort_site_dict = dict(sorted(site_dict.items()))
    # outfile="/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_development/crispr_off_target/20230831_demo_test/test.stat.xls"
    with open(outfile, "w") as w:
        header_list = ["Bulge_size", "Mismatches", "Sites", "DNAs"]
        w.write("\t".join(header_list) + "\n")
        for key, value in sort_site_dict.items():
            bulge_size = str(key[0])
            mismatch = str(key[1])
            sites = str(len(set(value)))
            dnas = str(len(value))
            new_line_list = [bulge_size, mismatch, sites, dnas]
            w.write("\t".join(new_line_list) + "\n")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="table统计")
    parser.add_argument("-i", "--infile", help="table表", required=True)
    parser.add_argument("-o", "--outfile", help="输出文件", required=True)
    args = parser.parse_args()
    stat_homo_region(args.infile, args.outfile)
