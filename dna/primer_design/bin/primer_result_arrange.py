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

def arranger_primer_result(infile: str, outfile: str) -> None:
    #infile = "/mnt/lustre/users/sanger-dev/sg-users/yuan.xu/majorbio_task/BSA/MJ20230620131_FX2023070500173_周玉峰_BSA/20230714_primer_design/primer_design_result/07.result/gprime.variant.result.xls"
    with open(infile, "r") as f, open(outfile, "w") as w:
        header_list = ["Chrom", "Pos", "Vtype", "Ref", "Alt", "Marker size(bp)", "Marker start(bp)", "Marker end(bp)", "FORWARD PRIMER (5'-3')",
                       "Tm(.C)", "GC(%)", "size", "REVERSE PRIMER (5'-3')", "Tm(.C)", "GC(%)", "size", "PRODUCT size (bp)", "start (bp)",
                       "end (bp)"]
        w.write("\t".join(header_list))
        w.write("\n")
        for line in f:
            if line.startswith("#"):
                continue
            line_list = line.strip().split("\t")
            if len(line_list) < 42:
                continue
            chrom, pos, total_number, vtype, ref, alt, mark_size, marker_start, marker_end, forward_primer_1, tm_f_1, gc_f_1, size_f_1, reverse_primer_1, tm_r_1, gc_r_1, size_r_1, product_size_1, product_start_1, product_end_1, forward_primer_2, tm_f_2, gc_f_2, size_f_2, reverse_primer_2, tm_r_2, gc_r_2, size_r_2, product_size_2, product_start_2, product_end_2, forward_primer_3, tm_f_3, gc_f_3, size_f_3, reverse_primer_3, tm_r_3, gc_r_3, size_r_3, product_size_3, product_start_3, product_end_3 = line_list
            new_line_list_1 = [chrom, pos, vtype, ref, alt, mark_size, marker_start, marker_end, forward_primer_1, tm_f_1, gc_f_1, size_f_1, reverse_primer_1, tm_r_1, gc_r_1, size_r_1, product_size_1, product_start_1, product_end_1]
            new_line_list_2 = [chrom, pos, vtype, ref, alt, mark_size, marker_start, marker_end, forward_primer_2, tm_f_2, gc_f_2, size_f_2, reverse_primer_2, tm_r_2, gc_r_2, size_r_2, product_size_2, product_start_2, product_end_2]
            new_line_list_3 = [chrom, pos, vtype, ref, alt, mark_size, marker_start, marker_end, forward_primer_3, tm_f_3, gc_f_3, size_f_3, reverse_primer_3, tm_r_3, gc_r_3, size_r_3, product_size_3, product_start_3, product_end_3]
            w.write("\t".join(new_line_list_1))
            w.write("\n")
            w.write("\t".join(new_line_list_2))
            w.write("\n")
            w.write("\t".join(new_line_list_3))
            w.write("\n")
            
            

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="整理primer3的结果")
    parse.add_argument("-i", "--infile", help="输入的primer3的结果", required=True)
    parse.add_argument("-o", "--outfile", help="输出整理好的结果文件", required=True)
    args = parse.parse_args()
    arranger_primer_result(args.infile, args.outfile)
