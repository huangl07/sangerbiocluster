# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20230307
# last_modified: 20230307


import argparse
from itertools import islice


def process_indel_anno_result(group_file: str, indel_anno_stat: str, indel_anno_stat_bsa: str) -> None:
    mix_info_dict = {}
    with open(group_file, "r", encoding='utf-8') as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                sample_name, mix_info, min_depth, max_depth = line_list
                mix_info_dict[mix_info] = sample_name
            else:
                raise Exception("group_file不为4列,所以报错")
    with open(indel_anno_stat, "r", encoding='utf-8') as f, open(indel_anno_stat_bsa, "w", encoding='utf-8') as m:
        for line in islice(f, 0, 1):
            line_list = line.strip().split("\t")
            m.write("\t".join(line_list))
            m.write("\n")
        for line in islice(f, 0, None):
            line_list = line.strip().split("\t")
            if line_list[0] in list(mix_info_dict.values()):
                m.write("\t".join(line_list))
                m.write("\n")

parser = argparse.ArgumentParser(description="生成snp_anno_stat.xls")
parser.add_argument("-g", "--group_file", required=True)
parser.add_argument("-s", "--indel_anno_stat", required=True)
parser.add_argument("-t", "--indel_anno_stat_bsa", required=True)
args = vars(parser.parse_args())
process_indel_anno_result(args["group_file"], args["indel_anno_stat"], args["indel_anno_stat_bsa"])
