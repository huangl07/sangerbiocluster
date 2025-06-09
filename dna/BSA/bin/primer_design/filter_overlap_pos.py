# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/05/24
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import vcfpy

def get_snp_list(infile: str, outfile: str, no_overlap_len: int) -> None:
    # 读取文件(支持.gz形式),最后一个会没有
    reader = vcfpy.Reader.from_path(infile)
    threshold_down = 0 
    threshold_up = 0
    overlap_count = 0
    with open(outfile, "w") as w:
        for record in reader:
            if threshold_down < record.POS < threshold_up:
                overlap_count += 1
                continue
            else: 
                if overlap_count == 1:
                    w.write("\t".join(line_list))
                    w.write("\n")
                overlap_count = 1
                threshold_down = record.POS - no_overlap_len
                threshold_up = record.POS + no_overlap_len
                line_list = [str(record.CHROM), str(record.POS), record.REF]
                line_list += [",".join([alt.value for alt in record.ALT])]
        if overlap_count == 1:
            w.write("\t".join(line_list))
            w.write("\n")
        else:
            pass
        

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="从vcf文件中获取snp_list文件")
    parse.add_argument("-i", "--infile", help="输入的vcf格式文件", required=True)
    parse.add_argument("-o", "--outfile", help="输出snp_pos.list", required=True)
    parse.add_argument("-l", "--no_overlap_len", type = int, help="snp不重叠区域的大小", required=True)
    args = parse.parse_args()
    get_snp_list(args.infile, args.outfile, args.no_overlap_len)
