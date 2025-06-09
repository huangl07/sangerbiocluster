# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20230307
# last_modified: 20230307

"""
通过qc.stat和group_list文件生成rawdata.xls和qc.xls
"""

import argparse
from itertools import islice


def process_qc_result(qc_stat: str, group_file: str, raw_data: str, clean_data: str) -> None:
    mix_info_dict = {}
    with open(group_file, "r", encoding='utf-8') as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                sample_name, mix_info, min_depth, max_depth = line_list
                mix_info_dict[mix_info] = sample_name
            else:
                raise Exception("group_file不为4列,所以报错")
    with open(qc_stat, "r", encoding='utf-8') as f, open(raw_data, "w", encoding='utf-8') as m, open(clean_data, "w", encoding='utf-8') as n:
        raw_header_list = ["Sample ID", "Raw Reads", "Raw Bases（bp）", "Raw GC（%）", "Raw Q30（%）"]
        clean_header_list = ["Sample ID", "Clean Reads", "Clean Bases（bp）", "Clean GC（%）", "Clean Q30（%）"]
        m.write("\t".join(raw_header_list))
        m.write("\n")
        n.write("\t".join(clean_header_list))
        n.write("\n")
        for line in islice(f, 1, None):
            line_list = line.strip().split("\t")
            if line_list[0] in list(mix_info_dict.values()):
                sample_id = line_list[0]
                raw_reads = line_list[1]
                raw_bases = line_list[2]
                raw_gc = round(float(line_list[5])*100, 2).__str__()
                raw_q30 = round(float(line_list[4])*100, 2).__str__()
                clean_reads = line_list[7]
                clean_bases = line_list[8]
                clean_gc = round(float(line_list[11])*100, 2).__str__()
                clean_q30 = round(float(line_list[10])*100, 2).__str__()
                raw_new_line = [sample_id, raw_reads, raw_bases, raw_gc, raw_q30]
                clean_new_line = [sample_id, clean_reads, clean_bases, clean_gc, clean_q30]
                m.write("\t".join(raw_new_line))
                m.write("\n")
                n.write("\t".join(clean_new_line))
                n.write("\n")

parser = argparse.ArgumentParser(description="通过qc.stat和group_list文件生成rawdata.xls和qc.xls")
parser.add_argument("-q", "--qc_stat", required=True)
parser.add_argument("-g", "--group_file", required=True)
parser.add_argument("-r", "--raw_data", required=True)
parser.add_argument("-c", "--clean_data", required=True)
args = vars(parser.parse_args())
process_qc_result(args["qc_stat"], args["group_file"], args["raw_data"], args["clean_data"])
