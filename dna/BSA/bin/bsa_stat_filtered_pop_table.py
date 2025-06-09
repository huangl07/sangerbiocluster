# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20220222
# last_modified: 20220228

"""
pop.table.filtered统计
"""

import argparse
import pandas as pd


def stat_group(group: pd.DataFrame) -> pd.DataFrame:
    """
    统计总的位点数量和有限的位点数量（snpeff注释结果为HIGH或者MODERATE）
    """
    number_total = group.shape[0]
    number_high = group[group["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[2] == "HIGH"].shape[0]
    number_moderate = group[group["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[2] == "MODERATE"].shape[0]
    effective_number = number_high + number_moderate
    new_df = pd.concat([pd.Series(number_total), pd.Series(effective_number)], axis=1)
    return new_df

def stat_gene(group: pd.DataFrame) -> pd.DataFrame:
    """
    统计pop.table.filtered总的位点数量和effective的位点数量（snpeff注释结果为HIGH或者MODERATE）
    """
    gene_total = group["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[3].drop_duplicates(keep="first")
    gene_total_number = gene_total.shape[0]
    gene_high_df = group[group["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[2] == "HIGH"]
    gene_moderate_df = group[group["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[2] == "MODERATE"]
    if not gene_high_df.empty:
        gene_high = gene_high_df["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[3].drop_duplicates(keep="first")
        gene_high_number = gene_high.shape[0]
    else:
        gene_high_number = 0
    if not gene_moderate_df.empty:
        gene_moderate = gene_moderate_df["ANN"].str.split(";", expand=True)[0].str.split("|", expand=True)[3].drop_duplicates(keep="first")
        gene_moderate_number = gene_moderate.shape[0]
    else:
        gene_moderate_number = 0
    gene_effective_number = gene_high_number + gene_moderate_number
    new_df = pd.concat([pd.Series(gene_total_number), pd.Series(gene_effective_number)], axis=1)
    return new_df

def stat_filtered_pop_table(infile: str, outfile: str , chrfile: str) -> None:
    """
    主函数
    """
    chr_list = pd.read_csv(chrfile, sep="\t", header=None)[0].values.tolist()
    df = pd.read_csv(infile, sep="\t")
    df = df.loc[df['CHROM'].isin(chr_list),]
    gene_df = df.groupby("CHROM").apply(stat_gene)
    gene_df.index = [index[0] for index in gene_df.index] # 多级索引降级，写在函数内会报错
    gene_df.reset_index(inplace=True) # 索引变第一列
    gene_df.rename({"index": "CHROM", 0:"Gene_Number", 1:"Effective_Gene"}, axis=1, inplace=True)
    snp_df = df[df["Vtype"] == "SNP"]
    indel_df = df[df["Vtype"] == "INDEL"]
    snp_stat = snp_df.groupby("CHROM").apply(stat_group)
    snp_stat.index = [index[0] for index in snp_stat.index] # 多级索引降级，写在函数内会报错
    snp_stat.reset_index(inplace=True) # 索引变第一列
    snp_stat.rename({"index": "CHROM", 0:"SNP_Number", 1:"Effective_SNP"}, axis=1, inplace=True)
    indel_stat = indel_df.groupby("CHROM").apply(stat_group)
    indel_stat.index = [index[0] for index in indel_stat.index]
    indel_stat.reset_index(inplace=True)
    indel_stat.rename({"index": "CHROM", 0:"INDEL_Number", 1:"Effective_INDEL"}, axis=1, inplace=True)
    merge_df = pd.merge(pd.DataFrame(snp_stat), pd.DataFrame(indel_stat), how="outer", on="CHROM", sort=True).\
        fillna(0).astype({"CHROM": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    merge_df = pd.merge(pd.DataFrame(merge_df), pd.DataFrame(gene_df), how="outer", on="CHROM", sort=True).\
        fillna(0).astype({"CHROM": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int,
                          "Gene_Number": int, "Effective_Gene": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    merge_df.to_csv(outfile, sep="\t", index=False)


if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="pop.table.filtered文件统计")
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-c", "--chrfile", required=True)
    args = vars(parser.parse_args())
    stat_filtered_pop_table(args["infile"], args["outfile"], args["chrfile"])