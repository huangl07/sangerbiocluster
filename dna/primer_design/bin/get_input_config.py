# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/05/24
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse

def get_input_config(infile: str, outfile: str, snp_list: str, fa_length: int) -> None:
    snp_dict = {}
    with open(snp_list, "r") as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                chrom, pos, ref, alt = line_list
                snp_dict[chrom + "-" + pos] = "[" + ref + "/" + alt + "]"
            else:
                raise Exception("输入文件不为4列")
    with open(infile, "r") as f, open(outfile, "w") as w:
        for line in f:
            if line.startswith(">"):
                chrom, bed = line.strip().split(":")
                chrom = chrom.replace(">", "")
                start_pos, end_pos = bed.strip().split("-")
                pos = int((int(start_pos) + int(end_pos))/2)
                config = chrom + "-" + str(pos)
                slice_pos = int(pos) - int(start_pos)
            else:
                seq = line.strip()
                if len(seq) < fa_length:
                    continue
                new_seq = seq[0:slice_pos] + snp_dict[config] + seq[slice_pos+1:None]
                new_line_list = [config, chrom, new_seq]
                w.write(",".join(new_line_list))
                w.write("\n")
                
        
if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="获取位点上下游的一段位置")
    parse.add_argument("-i", "--infile", help="输入fasta文件", required=True)
    parse.add_argument("-o", "--outfile", help="输出t_polymarker_input.csv", required=True)
    parse.add_argument("-s", "--snp_list", help="输入snp_list文件", required=True)
    parse.add_argument("-l", "--fa_length", help="输入fa长度", default=101)
    args = parse.parse_args()
    get_input_config(args.infile, args.outfile, args.snp_list, args.fa_length)
