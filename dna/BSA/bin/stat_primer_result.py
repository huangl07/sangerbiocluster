# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/07/06
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

"""
BSA引物设计结果
TODO:如果某个文件为空，还没考虑
"""

import argparse
import pandas as pd
import os

def stat_caps_result(caps_file: str) -> pd.DataFrame:
    """
    统计caps统计结果
    """
    caps_df = pd.read_csv(caps_file, sep="\t", header=0)
    caps_count = caps_df["Pos"].value_counts().shape[0]
    new_df = pd.concat([pd.Series(caps_count)], axis=1)
    new_df.rename({0:"caps_count"}, axis=1, inplace=True)
    return new_df

def stat_dcaps_result(dcaps_file: str) -> pd.DataFrame:
    dcaps_df = pd.read_csv(dcaps_file, sep="\t", header=0)
    dcaps_count = dcaps_df["Pos"].value_counts().shape[0]
    new_df = pd.concat([pd.Series(dcaps_count)], axis=1)
    new_df.rename({0:"dcaps_count"}, axis=1, inplace=True)
    return new_df

def stat_kasp_result(kasp_file: str) -> pd.DataFrame:
    """
    统计kasp的结果
    """
    kasp_df = pd.read_csv(kasp_file, sep="\t", header=0)
    kasp_count = kasp_df["position"].value_counts().shape[0]
    new_df = pd.concat([pd.Series(kasp_count)], axis=1)
    new_df.rename({0:"kasp_count"}, axis=1, inplace=True)
    return new_df

def stat_primer3_result(variant_snp_file: str, variant_indel_file: str) -> pd.DataFrame:
    """
    统计primer3结果
    """
    primer_snp_df = pd.read_csv(variant_snp_file, sep="\t", header=0)
    primer_indel_df = pd.read_csv(variant_indel_file, sep="\t", header=0)
    primer_snp_count = primer_snp_df["Pos"].value_counts().shape[0]
    primer_indel_count = primer_indel_df["Pos"].value_counts().shape[0]
    new_df = pd.concat([pd.Series(primer_snp_count+primer_indel_count)], axis=1)
    new_df.rename({0:"primer_count"}, axis=1, inplace=True)
    return new_df
  
def stat_one_method_result(method: str, caps_file: str, dcaps_file: str, kasp_file: str, variant_snp_file: str, 
                           variant_indel_file: str) -> pd.DataFrame:
    """
    统计一个方法的引物设计结果
    """
    caps_count_df = stat_caps_result(caps_file)
    dcaps_count_df = stat_dcaps_result(dcaps_file)
    kasp_count_df = stat_kasp_result(kasp_file)
    primer3_count_df = stat_primer3_result(variant_snp_file, variant_indel_file)
    method_count_df = pd.concat([caps_count_df, dcaps_count_df, kasp_count_df, primer3_count_df], axis=1)
    method_count_df.reset_index(inplace=True)
    method_count_df.rename({"index": "Method", "caps_count": "CAPS Count", "dcaps_count": "dCAPS Count",
                            "kasp_count": "KASP Count", "primer_count": "Sanger Count"}, axis=1, inplace=True)
    method_count_df["Method"][0] = method
    return method_count_df

def stat_primer_design_result(indir: str, outfile: str) -> None:
    """
    主函数
    """
    method_result_df = pd.DataFrame()
    for dir in os.listdir(indir):
        if os.path.isdir(os.path.join(indir, dir)):
            method = dir
            caps_file = os.path.join(indir, dir, method + ".caps.result.xls")
            dcaps_file = os.path.join(indir, dir, method + ".dcaps.result.xls")
            kasp_file = os.path.join(indir, dir, method + ".kasp.result.xls")
            variant_snp_file = os.path.join(indir, dir, method + ".sanger.snp.result.xls")
            variant_indel_file = os.path.join(indir, dir, method + ".sanger.indel.result.xls")
            one_method_result = stat_one_method_result(method, caps_file, dcaps_file, kasp_file, variant_snp_file, variant_indel_file)
            method_result_df = pd.concat([method_result_df, one_method_result], axis=0)
    method_result_df.to_csv(outfile, sep="\t", index=False)
            
if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="统计primer_design的结果")
    parser.add_argument("-i", "--indir", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    args = vars(parser.parse_args())
    stat_primer_design_result(args["indir"], args["outfile"])