# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20220222
# last_modified: 20220228

"""
处理deepBSA的结果
"""

import argparse
import pandas as pd

def process_deepBSA(infile: str, outfile: str) -> None:
    """
    主函数
    """
    df = pd.read_csv(infile, sep="\t", names=["X.chr", "pos", "value", "loess", "CI"])
    df[["pos"]] = df[["pos"]]/1000000
    df[["pos"]] = df[["pos"]].astype(int)
    df.to_csv(outfile, sep="\t", index=False)

if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="处理deepBSA的结果")
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    args = vars(parser.parse_args())
    process_deepBSA(args["infile"], args["outfile"])