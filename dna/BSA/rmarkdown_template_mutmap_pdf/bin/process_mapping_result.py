# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20230307
# last_modified: 20230307

"""
通过all.summary.stats和group_list文件生成align_stat.xls
"""

import argparse
from itertools import islice


def process_mapping_result(all_summary_stats: str, group_file: str, bsa_align_stat: str, bsa_coverage_stat: str) -> None:
    mix_info_dict = {}
    with open(group_file, "r", encoding='utf-8') as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                sample_name, mix_info, min_depth, max_depth = line_list
                mix_info_dict[mix_info] = sample_name
            else:
                raise Exception("group_file不为4列,所以报错")
    with open(all_summary_stats, "r", encoding='utf-8') as f, open(bsa_align_stat, "w", encoding='utf-8') as m, open(bsa_coverage_stat, "w", encoding='utf-8') as n:
        mapping_header_list = ["Sample ID", "Mapped Ratio（%）", "Properly Mapped（%）", "Duplication Ratio（%）"]
        coverage_header_list = ["Sample ID", "Average Insert Size", "Average Depth", "Coverage (>=1x)", "Coverage (>=4x)"]
        m.write("\t".join(mapping_header_list))
        m.write("\n")
        n.write("\t".join(coverage_header_list))
        n.write("\n")
        for line in islice(f, 1, None):
            line_list = line.strip().split("\t")
            if line_list[0] in list(mix_info_dict.values()):
                sample_id = line_list[0]
                mapping_ratio = line_list[1]
                proper_ratio = line_list[2]
                duplicate_ratio = line_list[3]
                insert_size = line_list[4]
                average_depth = line_list[5]
                coverage1 = line_list[7]
                coverage4 = line_list[8]
                mapping_new_line = [sample_id, mapping_ratio, proper_ratio, duplicate_ratio]
                coverage_new_line = [sample_id, insert_size, average_depth, coverage1, coverage4]
                m.write("\t".join(mapping_new_line))
                m.write("\n")
                n.write("\t".join(coverage_new_line))
                n.write("\n")

parser = argparse.ArgumentParser(description="通过all.summary.stats和group_list文件生成align_stat.xls和coverage_sample.xls")
parser.add_argument("-a", "--all_summary_stats", required=True)
parser.add_argument("-g", "--group_file", required=True)
parser.add_argument("-b", "--bsa_align_stat", required=True)
parser.add_argument("-c", "--bsa_coverage_stat", required=True)
args = vars(parser.parse_args())
process_mapping_result(args["all_summary_stats"], args["group_file"], args["bsa_align_stat"], args["bsa_coverage_stat"])
