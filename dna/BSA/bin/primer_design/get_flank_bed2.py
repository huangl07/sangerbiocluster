# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/05/29
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import os

def get_flank_bed(infile: str) -> None:
    if not os.path.isdir("flank_bed") and not os.path.exists("flank_bed"):
        os.makedirs("flank_bed")
    else:
        pass
    snp_id_dict = {}
    with open(infile, "r") as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                snp_id, chrom, pos, strand = line_list
                start_pos, end_pos = pos.split("-")
                snp_id_dict.setdefault(snp_id, []).append((chrom, start_pos, end_pos))
            else:
                raise Exception("输入文件不为4列")
    for key,value in snp_id_dict.items():
        outfile = os.path.join("flank_bed", key + ".txt")
        with open(outfile, "w") as w:
            for chrom, start_pos, end_pos in value:
                new_line_list= [chrom, start_pos, end_pos]
                w.write("\t".join(new_line_list))
                w.write("\n")

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="获取位点上下游的一段位置")
    parse.add_argument("-i", "--infile", help="输入的flank_ranges文件", required=True)
    args = parse.parse_args()
    get_flank_bed(args.infile)
