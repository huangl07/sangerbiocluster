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
    duplicate_df = group.drop_duplicates(subset="POS",keep="first")
    number_total = duplicate_df.shape[0]
    number_high = duplicate_df[duplicate_df["PutativeImpact"] == "HIGH"].shape[0]
    number_moderate = duplicate_df[duplicate_df["PutativeImpact"] == "MODERATE"].shape[0]
    effective_number = number_high + number_moderate
    new_df = pd.concat([pd.Series(number_total), pd.Series(effective_number)], axis=1)
    return new_df

def stat_gene(group: pd.DataFrame) -> pd.DataFrame:
    """
    统计pop.table.filtered总的位点数量和effective的位点数量（snpeff注释结果为HIGH或者MODERATE）
    """
    gene_total_number = group.shape[0]
    new_df = pd.concat([pd.Series(gene_total_number)], axis=1)
    return new_df

def stat_annotation(group: pd.DataFrame) -> pd.DataFrame:
    nrid_number = group[group["NRID"] != "--"].shape[0]
    Uniid_number = group[group["UniID"] != "--"].shape[0]
    koid_number = group[group["KoID"] != "--"].shape[0]
    goterm_number = group[group["GOTERM"] != "--"].shape[0]
    eggnog_number = group[group["EGGNOG"] != "--"].shape[0]
    pfamid_number = group[group["PfamID"] != "--"].shape[0]
    interproaccession_number = group[group["InterProAccession"] != "--"].shape[0]
    new_df = pd.concat([pd.Series(nrid_number), pd.Series(Uniid_number),pd.Series(koid_number),pd.Series(goterm_number),
                        pd.Series(eggnog_number),pd.Series(pfamid_number),pd.Series(interproaccession_number)], axis=1)
    return new_df

def stat_region_pop_table(infile: str, abstract_gene: str, anno_summary: str, outfile: str) -> None:
    """
    主函数
    """
    df = pd.read_csv(infile, sep="\t")
    genes_df = pd.read_csv(abstract_gene, sep="\t", dtype="object")
    anno_summary_df = pd.read_csv(anno_summary, sep="\t")
    new_anno_summary_df = pd.concat([anno_summary_df, anno_summary_df["GeneID"].str.split(":", expand=True)],axis=1).drop("GeneID", axis=1)
    new_anno_summary_df.rename(columns = {0: "gene_id", 1: "transcript_id", 2: "chr", 3: "gene_start", 4: "gene_end"}, inplace=True)
    merge_region_gene_df = pd.merge(genes_df, new_anno_summary_df, how="inner")
    gene_df = genes_df.groupby("REGION").apply(stat_gene)
    if genes_df.shape[0] == 1:
        gene_df.index = [genes_df["REGION"][0]]
    else:
        gene_df.index = [index[0] for index in gene_df.index] # 多级索引降级，写在函数内会报错
    gene_df.reset_index(inplace=True) # 索引变第一列
    gene_df.rename({"index": "REGION", 0:"Gene_Number", 1:"Effective_Gene"}, axis=1, inplace=True)
    snp_df = df[df["Vtype"] == "SNP"]
    indel_df = df[df["Vtype"] == "INDEL"]
    snp_stat = snp_df.groupby("REGION").apply(stat_group)
    snp_stat.index = [index[0] for index in snp_stat.index] # 多级索引降级，写在函数内会报错
    snp_stat.reset_index(inplace=True) # 索引变第一列
    snp_stat.rename({"index": "REGION", 0:"SNP_Number", 1:"Effective_SNP"}, axis=1, inplace=True)
<<<<<<< HEAD:Tumor/nextflow/bin/bsa_stat_region.py
    indel_stat = indel_df.groupby("REGION").apply(stat_group)
    indel_stat.index = [index[0] for index in indel_stat.index]
    indel_stat.reset_index(inplace=True)
    indel_stat.rename({"index": "REGION", 0:"INDEL_Number", 1:"Effective_INDEL"}, axis=1, inplace=True)
    merge_df = pd.merge(pd.DataFrame(snp_stat), pd.DataFrame(indel_stat), how="outer", on="REGION", sort=True).\
        fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    merge_df = pd.merge(pd.DataFrame(merge_df), pd.DataFrame(gene_df), how="outer", on="REGION", sort=True).\
        fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int,
                          "Gene_Number": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    Annotation_df = merge_region_gene_df.groupby("REGION").apply(stat_annotation)
    Annotation_df.index = [index[0] for index in Annotation_df.index] # 多级索引降级，写在函数内会报错
    Annotation_df.reset_index(inplace=True) # 索引变第一列
    Annotation_df.rename({"index": "REGION", 0:"NRID", 1:"UniID", 2:"KoID", 3:"GOTERM", 4:"EGGNOG",
                          5:"PfamID", 6:"InterProAccession"}, axis=1, inplace=True)
    merge_df = pd.merge(pd.DataFrame(merge_df), pd.DataFrame(Annotation_df), how="outer", on="REGION", sort=True).\
        fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int,
                          "Gene_Number": int, "NRID": int, "UniID": int, "KoID": int, "GOTERM": int,
                          "EGGNOG": int, "PfamID": int, "InterProAccession": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
=======
    if not indel_df.empty:
        indel_stat = indel_df.groupby("REGION").apply(stat_group)
        indel_stat.index = [index[0] for index in indel_stat.index]
        indel_stat.reset_index(inplace=True)
        indel_stat.rename({"index": "REGION", 0:"INDEL_Number", 1:"Effective_INDEL"}, axis=1, inplace=True)
        merge_df = pd.merge(pd.DataFrame(snp_stat), pd.DataFrame(indel_stat), how="outer", on="REGION", sort=True).\
            fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    else:
        merge_df = snp_stat.copy(deep=True)
        merge_df['INDEL_Number'] = 0
        merge_df['Effective_INDEL'] = 0
    if not gene_df.empty:
        merge_df = pd.merge(pd.DataFrame(merge_df), pd.DataFrame(gene_df), how="outer", on="REGION", sort=True).\
            fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int,
                              "Gene_Number": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
        Annotation_df = merge_region_gene_df.groupby("REGION").apply(stat_annotation)
        if merge_region_gene_df.shape[0] == 1:
            Annotation_df.index = [merge_region_gene_df["REGION"][0]]
        else:
            Annotation_df.index = [index[0] for index in Annotation_df.index] # 多级索引降级，写在函数内会报错
        Annotation_df.reset_index(inplace=True) # 索引变第一列
        Annotation_df.rename({"index": "REGION", 0:"NRID", 1:"UniID", 2:"KoID", 3:"GOTERM", 4:"EGGNOG",
                              5:"PfamID"}, axis=1, inplace=True)
        merge_df = pd.merge(pd.DataFrame(merge_df), pd.DataFrame(Annotation_df), how="outer", on="REGION", sort=True).\
            fillna(0).astype({"REGION": str, "SNP_Number": int, "Effective_SNP": int, "INDEL_Number": int, "Effective_INDEL": int,
                              "Gene_Number": int, "NRID": int, "UniID": int, "KoID": int, "GOTERM": int,
                              "EGGNOG": int, "PfamID": int}) # pd.merge会把int转变为float,同时缺失值变成NaN
    else:
        merge_df["Gene_Number"] = 0
        merge_df["NRID"] = 0
        merge_df["UniID"] = 0
        merge_df["KoID"] = 0
        merge_df["GOTERM"] = 0
        merge_df["EGGNOG"] = 0
        merge_df["PfamID"] = 0
>>>>>>> BSA_v1:BSA/bin/bsa_stat_region.py
    merge_df.to_csv(outfile, sep="\t", index=False)


if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="region区域信息统计文件统计")
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-a", "--abstract_gene", required=True)
    parser.add_argument("-b", "--anno_summary", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    args = vars(parser.parse_args())
    stat_region_pop_table(args["infile"], args["abstract_gene"], args["anno_summary"], args["outfile"])