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

def stat_primer3_result(snp_variant_file: str, indel_variant_file: str) -> pd.DataFrame:
    """
    统计primer3结果
    """
    snp_primer_df = pd.read_csv(snp_variant_file, sep="\t", header=0)
    indel_primer_df = pd.read_csv(indel_variant_file, sep="\t", header=0)
    snp_primer_count = snp_primer_df["Pos"].value_counts().shape[0]
    indel_primer_count = indel_primer_df["Pos"].value_counts().shape[0]
    new_df = pd.concat([pd.Series(snp_primer_count+indel_primer_count)], axis=1)
    new_df.rename({0:"primer_count"}, axis=1, inplace=True)
    return new_df
  
def stat_one_method_result(method: str, caps_file: str, dcaps_file: str, kasp_file: str, snp_variant_file: str, indel_variant_file: str) -> pd.DataFrame:
    """
    统计一个方法的引物设计结果
    caps_count_df = stat_caps_result(caps_file)
    dcaps_count_df = stat_dcaps_result(dcaps_file)
    kasp_count_df = stat_kasp_result(kasp_file)
    primer3_count_df = stat_primer3_result(snp_variant_file, indel_variant_file)
    method_count_df = pd.concat([caps_count_df, dcaps_count_df, kasp_count_df, primer3_count_df], axis=1)
    method_count_df.reset_index(inplace=True)
    method_count_df.rename({"index": "Method", "caps_count": "CAPS Count", "dcaps_count": "dCAPS Count",
                            "kasp_count": "KASP Count", "primer_count": "Sanger Count"}, axis=1, inplace=True)
    method_count_df["Method"][0] = method
    return method_count_df
    """
    method_count_df = pd.DataFrame({"Method": [method]})
    
    if os.path.exists(caps_file):
        caps_count_df = stat_caps_result(caps_file)
    else:
        caps_count_df = pd.DataFrame({"caps_count": [0]})
    method_count_df = pd.concat([method_count_df, caps_count_df], axis=1)
    
    if os.path.exists(dcaps_file):
        dcaps_count_df = stat_dcaps_result(dcaps_file)
    else:
        dcaps_count_df = pd.DataFrame({"dcaps_count": [0]})
    method_count_df = pd.concat([method_count_df, dcaps_count_df], axis=1)
    
    if os.path.exists(kasp_file):
        kasp_count_df = stat_kasp_result(kasp_file)
    else:
        kasp_count_df = pd.DataFrame({"kasp_count": [0]})
    method_count_df = pd.concat([method_count_df, kasp_count_df], axis=1)
    
    if os.path.exists(snp_variant_file) and os.path.exists(indel_variant_file):
        primer3_count_df = stat_primer3_result(snp_variant_file, indel_variant_file)
    else:
        primer3_count_df = pd.DataFrame({"primer_count": [0]})
    method_count_df = pd.concat([method_count_df, primer3_count_df], axis=1)
    
    method_count_df.rename({"caps_count": "CAPS Count", "dcaps_count": "dCAPS Count", "kasp_count": "KASP Count", "primer_count": "Sanger Count"}, axis=1, inplace=True)
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
            snp_variant_file = os.path.join(indir, dir, method + ".sanger.snp.result.xls")
            indel_variant_file = os.path.join(indir, dir, method + ".sanger.snp.result.xls")
            one_method_result = stat_one_method_result(method, caps_file, dcaps_file, kasp_file, snp_variant_file, indel_variant_file)
            method_result_df = pd.concat([method_result_df, one_method_result], axis=0)
    method_result_df.to_csv(outfile, sep="\t", index=False)
            
if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="统计primer_design的结果")
    parser.add_argument("-i", "--indir", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    args = vars(parser.parse_args())
    stat_primer_design_result(args["indir"], args["outfile"])
