# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/08/25
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import pandas as pd


def separate_table(vcf_table, outfile):
    #vcf_table = "/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_project/shouhou/xiaofeidi_vcf2table/finalvcf/BCO9_tissue_vs_NA_somatic_final_anno.table.xls"
    vcf_table_df = pd.read_csv(vcf_table, sep="\t", encoding="utf-8")
    vcf_table_df["ANN"] = vcf_table_df["ANN"].str.split(";")
    vcf_table_df = vcf_table_df.explode("ANN", ignore_index=True)
    # 兼容vep注释
    if len(vcf_table_df["ANN"][0].split("|"))>16:
        anno_df = vcf_table_df["ANN"].str.split("|", expand=True).rename(columns={0: "Allele",
                                                                             1: "Consequence",
                                                                             2: "IMPACT",
                                                                             3: "SYMBOL",
                                                                             4: "GeneName",
                                                                             5: "FeatureType",
                                                                             6: "FeatureID",
                                                                             7: "BioType",
                                                                             8: "EXON",
                                                                             9: "INTRON",
                                                                             10: "HVGS.c",
                                                                             11: "HVGS.p",
                                                                             12: "cDNA_position",
                                                                             13: "CDS_position",
                                                                             14: "Protein_position",
                                                                             15: "Amino_acids",
                                                                             16: "Codons",
                                                                             17: "Existing_variation",
                                                                             18: "DISTANCE",
                                                                             19: "STRAND",
                                                                             20: "FLAGS",
                                                                             21: "SYMBOL_SOURCE",
                                                                             22: "HGNC_ID",
                                                                             23: "Source",
                                                                             24: "HGVS_offset",
                                                                             25: "gtf"
                                                                             })
    else:
        anno_df = vcf_table_df["ANN"].str.split("|", expand=True).rename(columns={0: "Allele",
                                                                             1: "Annotation",
                                                                             2:"PutativeImpact",
                                                                             3:"GeneName",
                                                                             4:"GeneID",
                                                                             5:"FeatureType",
                                                                             6:"FeatureID",
                                                                             7:"TranscriptBioType",
                                                                             8:"Rank",
                                                                             9:"HGVS.c",
                                                                             10:"HVGS.p",
                                                                             11:"cDNA_position_len",
                                                                             12:"CDS_position_len",
                                                                             13:"Protein_position_len",
                                                                             14:"Distance",
                                                                             15:"ERRORS"})
    variant_df = vcf_table_df.drop("ANN", axis=1)
    merge_df = pd.concat([variant_df, anno_df], axis=1)
    #outfile = "/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_project/shouhou/xiaofeidi_diff_variant/差异位点提取/BC109_T_vs_BC46_T.variant.diff.anno.xls"
    merge_df.to_csv(outfile, index=False, header=True, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="table整理")
    parser.add_argument("-i", "--vcf_table", help="vcf转成的table表", required=True)
    parser.add_argument("-o", "--outfile", help="输出文件", required=True)
    args = parser.parse_args()
    separate_table(args.vcf_table, args.outfile)
