# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20230307
# last_modified: 20230307


import argparse
from itertools import islice


def process_snp_result(group_file: str, snp_stat: str, snp_stat_bsa: str) -> None:
    mix_info_dict = {}
    with open(group_file, "r", encoding='utf-8') as f:
        for line in f:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                sample_name, mix_info, min_depth, max_depth = line_list
                mix_info_dict[mix_info] = sample_name
            else:
                raise Exception("group_file不为4列,所以报错")
    with open(snp_stat, "r", encoding='utf-8') as f, open(snp_stat_bsa, "w", encoding='utf-8') as m:
        snp_stat_bsa_header_list = ["Sample ID", "SNP Number", "Transition", "Transvertion", "Ts/Tv", "Heterozygosity Number", "Homozygosity Number"]
        m.write("\t".join(snp_stat_bsa_header_list))
        m.write("\n")
        for line in islice(f, 1, None):
            line_list = line.strip().split("\t")
            if line_list[0] in list(mix_info_dict.values()):
                sample_id = line_list[0]
                num_homo_ref = line_list[3]
                num_homo_alt = line_list[4]
                num_homo = int(num_homo_ref) + int(num_homo_alt)
                num_het = line_list[5]
                num_snv = line_list[6]
                num_ts = line_list[10]
                num_tv = line_list[11]
                ts_tv = line_list[12]
                new_line = [sample_id, num_snv, num_ts, num_tv, ts_tv, num_het, str(num_homo)]
                m.write("\t".join(new_line))
                m.write("\n")

parser = argparse.ArgumentParser(description="生成snp_stat.xls")
parser.add_argument("-g", "--group_file", required=True)
parser.add_argument("-s", "--snp_stat", required=True)
parser.add_argument("-b", "--snp_stat_bsa", required=True)
args = vars(parser.parse_args())
process_snp_result(args["group_file"], args["snp_stat"], args["snp_stat_bsa"])
