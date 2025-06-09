# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yiwei.tang
# first_modified: 20231226
# last_modified: 20230406

import argparse


def abstract_mix_info_from_group_file(
    group_file: str, mix_info: str, pop_info: str
) -> None:
    """
    主函数：从group_file文件中获取所有混池列表
    """
    with open(group_file, "r", encoding="utf-8") as f:
        with open(mix_info, "w", encoding="utf-8") as w, open(
            pop_info, "w", encoding="utf-8"
        ) as w2:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                line_list = line.strip().split("\t")
                if len(line_list) == 4:
                    sample_name, mix_info, _, _ = line.strip().split("\t")
                    w2.write(sample_name)
                    w2.write("\n")
                    if mix_info.startswith("B"):
                        w.write(sample_name)
                        w.write("\n")
                    else:
                        continue
                else:
                    raise Exception("group文件不为4列")


parser = argparse.ArgumentParser(description="从group_file文件中获取所有混池列表")
parser.add_argument("-i", "--group_file", required=True)
parser.add_argument("-o", "--mix_info", required=True)
parser.add_argument("-p", "--pop_info", required=True)
args = vars(parser.parse_args())
abstract_mix_info_from_group_file(
    args["group_file"], args["mix_info"], args["pop_info"]
)
