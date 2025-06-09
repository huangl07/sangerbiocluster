# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/08/23
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
from collections import OrderedDict
import pandas as pd
from typing import Dict, List, Mapping


def sample_list_extend(sample_list: List[str]) -> List[str]:
    """
    功能：扩展sample_list表格，例如：[B1, B2] => [B1.GT, B1.AD, B1.DP, B2.GT, B2.AD, B2.DP]
        ：param sample_list 1个包含样本名的列表
        ：return List[str]
    """
    new_sample_list = []
    for sample_name in sample_list:
        sample_name_GT = sample_name + '.GT'
        sample_name_AD = sample_name + '.AD'
        sample_name_DP = sample_name + '.DP'
        new_sample_list.append(sample_name_GT)
        new_sample_list.append(sample_name_AD)
        new_sample_list.append(sample_name_DP)
    return new_sample_list


def sample_dict_extend(sample_dict: Mapping[str, str]) -> Dict[str, str]:
    """
    功能：把sample_dict扩充，例如：{"SRR5516308": "B1", "SRR5516308": "B2"} => \
        {"SRR5516308_GT": "B1", "SRR551630_AD": "B1", "SRR551630_DP": "B1", \
        "SRR5516308_GT": "B2", "SRR5516308_AD": "B2", "SRR5516308_DP": "B2"}
        ：param sample_dict 一个需要修改的的字典
        : return Dict[str, str]
    """
    new_sample_dict = {}
    for key, value in sample_dict.items():
        new_sample_dict[key + '.GT'] = value + '.GT'
        new_sample_dict[key + '.AD'] = value + '.AD'
        new_sample_dict[key + '.DP'] = value + '.DP'
    return new_sample_dict


def rename_df(df: pd.DataFrame, sample_detail_dict: Mapping[str, str]) -> pd.DataFrame:
    """
    功能：重命名列名
        :param df
        :param 要重命名的字典
        :return df
    """
    return df.rename(columns=sample_detail_dict)


def filter_polyallelic(df: pd.DataFrame) -> pd.DataFrame:
    """
    功能:多等位基因位点过滤，df的Alt列
        :param df
        :returm df
    """
    return df[~ df["Alt"].str.contains(",")]


def filter_by_dp(df: pd.DataFrame, lowdepth : int) -> pd.DataFrame:
    """
    功能：根据深度阈值过滤位点(空的df.all()为True)
        ：param df
        : param lowdepth 最低深度阈值
        ：return df
    """
    filtered = df[(df.filter(regex = '.*\.DP', axis=1) > lowdepth).all(axis=1)]
    result = filtered
    return result

def filter_by_gt(df: pd.DataFrame) -> pd.DataFrame:
    """
    功能: 过滤基因型一样的位点
        ：param df
        ：return df
    """
    filtered = df[(df.filter(regex = r'.*\.GT', axis=1).apply(lambda x: len(set(x))==len(x),axis=1))]
    result = filtered
    return result

def main(table_file: str, group_file: str, lowdepth: int, outfile: str) -> None:
    """
    主函数
    01.抽取所需要的列（亲本和混池）
    02.修改列名
    03.常规过滤，过滤基因型为None/None(python结果)和./.(gatk结果)的位点
    04.多等位基因过滤
    """
    sample_list = []
    sample_dict = OrderedDict()
    with open(group_file, "r") as f:
        for line in f:
            sample_id, sample_info = line.strip().split("\t")[:2]
            if sample_id != "-":
                sample_list.append(sample_id)
                sample_dict[sample_id] = sample_info
            else:
                pass
    df = pd.read_csv(table_file, sep="\t")  # 读取vcf_table文件

    sample_detail_list = sample_list_extend(sample_list)
    sample_detail_dict = sample_dict_extend(sample_dict)
    abstract_list = ["CHROM", "POS", "Ref", "Alt", "Vtype"] + sample_detail_list + ["ANN"]
    abstract_df = df[abstract_list]
    new_abstract_df1 = abstract_df[(abstract_df[sample_detail_list] != "None/None").all(axis=1)]
    new_abstract_df2 = new_abstract_df1[(new_abstract_df1[sample_detail_list] != "./.").all(axis=1)]
    new_abstract_df3 = filter_polyallelic(rename_df(new_abstract_df2, sample_detail_dict))
    new_abstract_df = new_abstract_df3
    """
    根据传入的深度阈值过滤
    """
    filtered_df2 = filter_by_gt(filter_by_dp(new_abstract_df, lowdepth))
    result_df = filtered_df2
    result_df.to_csv(outfile, sep="\t", index=False)

if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="pop.table根据分组文件过滤")
    parser.add_argument("-i", "--table_file", required=True)
    parser.add_argument("-g", "--group_file", required=True)
    parser.add_argument("-l", "--lowdepth", type=int, required=True, default=10)
    parser.add_argument("-o", "--outfile", required=True)
    args = vars(parser.parse_args())

    main(args["table_file"], args["group_file"], args["lowdepth"], args["outfile"])

