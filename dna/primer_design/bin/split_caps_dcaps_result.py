# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/07/19
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import csv
from collections import namedtuple
from itertools import islice 

def split_caps_dcaps_result(infile: str, caps_result: str, dcaps_result: str) -> None:
    #infile = "/mnt/lustre/users/sanger-dev/sg-users/yuan.xu/majorbio_task/BSA/MJ20230620131_FX2023070500173_周玉峰_BSA/20230714_primer_design/primer_design_result/07.result/gprime.caps_primers_result.xls"
    with open(infile, "r") as f, open(caps_result, "w") as m, open(dcaps_result, "w") as n:
        header_list = ["Chrom", "Pos", "Ptype", "Enzyme", "Bases", "Product_size", "TM_left", "GCcontent_left", "Primer_seq_left",
                       "TM_right", "GCcontent_right", "Primer_seq_right"]
        m.write("\t".join(header_list))
        m.write("\n")
        n.write("\t".join(header_list))
        n.write("\n")
        for line in islice(f, 0, 1):
            pass
        for line in islice(f, 0, None):
            line_list = line.strip().split("\t")
            chrom, pos, ptype, enzyme, bases, product_size, tm_left, gccontent_left, primer_seq_left, tm_right, gccontent_right, primer_seq_right = line_list
            if ptype == "CAPS":
                m.write(line)
            elif ptype == "dCAPS":
                n.write(line)
            

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="整理primer3的结果")
    parse.add_argument("-i", "--infile", help="输入的vcf格式文件", required=True)
    parse.add_argument("-c", "--caps_result", help="输出caps_result结果", required=True)
    parse.add_argument("-d", "--dcaps_result", help="输出dcaps_result结果", required=True)
    args = parse.parse_args()
    split_caps_dcaps_result(args.infile, args.caps_result, args.dcaps_result)
