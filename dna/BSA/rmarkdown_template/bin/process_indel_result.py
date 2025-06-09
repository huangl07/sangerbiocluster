# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20230307
# last_modified: 20230307


import argparse
from itertools import islice


def process_indel_result(group_file: str, indel_stat: str, indel_stat_bsa: str) -> None:
    mix_info_dict = {}
    with open(group_file, "r", encoding='utf-8') as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                sample_name, mix_info, min_depth, max_depth = line_list
                mix_info_dict[mix_info] = sample_name
            else:
                raise Exception("group_file不为4列,所以报错")
    with open(indel_stat, "r", encoding='utf-8') as f, open(indel_stat_bsa, "w", encoding='utf-8') as m:
        indel_stat_bsa_header_list = ["Sample ID", "Insert Number", "Delete Number", "Heterozygosity Number", "Homozygosity Number"]
        m.write("\t".join(indel_stat_bsa_header_list))
        m.write("\n")
        for line in islice(f, 1, None):
            line_list = line.strip().split("\t")
            if line_list[0] in list(mix_info_dict.values()):
                sample_id = line_list[0]
                num_insert = line_list[1]
                num_delete = line_list[2]
                num_het = line_list[3]
                num_homo_ref = line_list[4]
                new_line = [sample_id, num_insert, num_delete, num_het, num_homo_ref]
                m.write("\t".join(new_line))
                m.write("\n")

parser = argparse.ArgumentParser(description="生成indel_stat.xls")
parser.add_argument("-g", "--group_file", required=True)
parser.add_argument("-s", "--indel_stat", required=True)
parser.add_argument("-b", "--indel_stat_bsa", required=True)
args = vars(parser.parse_args())
process_indel_result(args["group_file"], args["indel_stat"],  args["indel_stat_bsa"])
