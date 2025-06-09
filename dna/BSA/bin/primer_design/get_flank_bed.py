# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/05/24
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse

def get_flank_bed(infile: str, outfile: str, flank_len: int) -> None:
    with open(infile, "r") as f, open(outfile, "w") as w:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                chrom, pos, ref, alt = line_list
                new_line_list = [chrom, str(max(1, int(pos) - flank_len - 1)), str(int(pos) + flank_len)]
                w.write("\t".join(new_line_list))
                w.write("\n")
            else:
                raise Exception("输入文件不为4列")

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="获取位点上下游的一段位置")
    parse.add_argument("-i", "--infile", help="输入的table文件", required=True)
    parse.add_argument("-o", "--outfile", help="输出bed文件", required=True)
    parse.add_argument("-l", "--flank_len", type = int, help="上下游区域大小", required=True)
    args = parse.parse_args()
    get_flank_bed(args.infile, args.outfile, args.flank_len)
