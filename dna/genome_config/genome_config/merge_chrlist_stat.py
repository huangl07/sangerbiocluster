# -*- coding: utf-8 -*-
# @Last-edit Time 2022/12/7
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com

import pandas as pd
import argparse

parse = argparse.ArgumentParser(description="汇总chromosomelist")
parse.add_argument("-c", "--clist", help="输入chromosome_list", required=True)
parse.add_argument("-t", "--table", help="输入seqkit生成的序列统计", required=True)
parse.add_argument("-o", "--output", help="输出结果文件", required=True)
args = vars(parse.parse_args())

clist = pd.read_csv(args["clist"], sep="\t", header=None).iloc[:, [0, 1, 2, 2]]
clist.columns = ["initial_name", "assembly_level", "output_name", "cache_name"]
table = pd.read_csv(args["table"], sep="\t", header=0).iloc[:, [0, 1, 2]]
table.columns = ["output_name", "length", "GC"]
out = pd.merge(clist, table)
out.to_csv(args["output"], index=False, sep="\t")
