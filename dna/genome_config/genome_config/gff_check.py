# -*- coding: utf-8 -*-
# @Last-edit Time 2022/12/1
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com

import re
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="gff_check.py --gtf gff --output outputdir")
parser.add_argument("--gtf", type=str, default=None)
parser.add_argument("--output", type=str, default=".")
args = parser.parse_args()


def gene_filter(temp):
    gene_id = ""
    try:
        gene_id = temp.split('gene_id "')[1].split('"')[0]
    except:
        pass
    if gene_id:
        return gene_id
    return ""


def trans_filter(temp):
    trans_id = ""
    try:
        trans_id = temp.split('transcript_id "')[1].split('"')[0]
    except:
        pass
    if trans_id:
        return trans_id
    return ""


gtf_df = pd.read_csv(args.gtf, sep="\t", header=None, comment="#")
df_colname = [
    "chrom", "source", "type", "start", "end", "score", "straight", "frag",
    "attr"
]
gtf_df.columns = df_colname
gtf_df.insert(9, "gene_id", gtf_df.attr.apply(gene_filter))
gtf_df.insert(10, "trans_id", gtf_df.attr.apply(trans_filter))
pos_df = gtf_df.query(
    'gene_id != ""').loc[:, ["gene_id", "trans_id", "chrom", "start", "end"]]
pos_uniq_df = pos_df.groupby(["gene_id", "trans_id", "chrom"]).agg({
    "start": min,
    "end": max
})
pos_uniq_df.reset_index(inplace=True)
pos_uniq_df.columns = ["gene_id", "transcript_id", "chr", "start", "end"]
pos_uniq_df.to_csv(os.path.join(args.output, "gene_pos.txt"),
                   sep="\t",
                   header=True,
                   index=0)
pos_no_df = gtf_df.query('gene_id == ""').loc[:, df_colname]
pos_no_df.to_csv(os.path.join(args.output, "gene_pos_no.txt"),
                 sep="\t",
                 header=0,
                 index=0)
