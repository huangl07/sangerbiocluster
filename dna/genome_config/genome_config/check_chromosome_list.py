# coding: utf-8
"""检查chromosome_list
"""
import argparse
import pandas
import numpy

parse = argparse.ArgumentParser(description="生成chromosome_list")
parse.add_argument("-i", "--input", help="输入chromosome_list", required=True)
parse.add_argument("-f", "--fai", help="输入ref.fa.fai文件", required=True)
parse.add_argument("-o", "--output", help="输出结果", required=True)
args = vars(parse.parse_args())

df = pandas.read_csv(args["input"], sep="\t", header=None)
if df.shape[1] >= 3:
    out_df = df
else:
    df.columns = ["input", "level"]
    # 处理染色体部分
    chrdf = df.iloc[numpy.flatnonzero(df["level"] == "Chromosome"), ].copy()
    chrInt = chrdf.iloc[numpy.flatnonzero(
        chrdf.input.apply(lambda x: isinstance(x, int) or x.isdigit())
    ), ].copy()
    # 如果原始染色体名称包含整数
    if chrInt.shape[0] > 0:
        chrInt.insert(2, "intInput", chrInt.input.apply(int))
        chrInt.sort_values(by="intInput", inplace=True)
        chrInt.reset_index(drop=True, inplace=True)
        chrInt.insert(3, "outIndex", chrInt.index + 1)
        chrInt.insert(4, "out",
                      chrInt.outIndex.apply(lambda x: "chr" + str(x)))
        chrInt.drop(columns=['intInput', 'outIndex'], inplace=True)
        chrStr = chrdf.iloc[numpy.flatnonzero(
            chrdf.input.apply(lambda x: not isinstance(x, int) and not x.
                              isdigit())), ].copy()
        if chrStr.shape[0] > 0:
            chrStr.insert(
                2, "out",
                chrdf.input.apply(lambda x: "chr" + str(x).capitalize()))
        chrdf = pandas.concat([chrInt, chrStr])
    else:
        chrdf.insert(2, "outIndex", chrdf.index + 1)
        chrdf.insert(3, "out", chrdf.outIndex.apply(lambda x: "chr" + str(x)))
        chrdf.drop(columns=['outIndex'], inplace=True)
    sca = df.iloc[numpy.flatnonzero(df["level"] != "Chromosome"), ].copy()
    faindex = pandas.read_csv(args["fai"], header=None, sep="\t", index_col=0)
    sca.insert(2, "len", sca.input.apply(lambda x: faindex.iloc[:, 0][x]))
    sca.sort_values(by="len", ascending=False, inplace=True)
    sca.reset_index(drop=True, inplace=True)
    sca.insert(3, "outIndex", sca.index + 1)
    sca.insert(4, "out", sca.outIndex.apply(lambda x: "sca" + str(x)))
    sca.drop(columns=['len', 'outIndex'], inplace=True)
    out_df = pandas.concat([chrdf, sca])

out_df.columns = ["input", "level", "output"]
chrdf = out_df.iloc[numpy.flatnonzero(
    out_df["level"] == "Chromosome"), ].copy()
chrSort = chrdf.output.apply(lambda x: int(x[3:])
                             if x[3:].isdigit() else 99999)
chrdf.insert(3, "sort", chrSort)
chrdf.sort_values(by=["sort", "output"], inplace=True)
chrdf.reset_index(drop=True, inplace=True)
chrdf.drop(columns=["sort"], inplace=True)
sca = out_df.iloc[numpy.flatnonzero(out_df["level"] != "Chromosome"), ].copy()
scaSort = sca.output.apply(lambda x: int(x[3:]) if x[3:].isdigit() else 99999)
sca.insert(3, "sort", scaSort)
sca.sort_values(by=["sort", "output"], inplace=True)
sca.reset_index(drop=True, inplace=True)
sca.drop(columns=["sort"], inplace=True)
out_df = pandas.concat([chrdf, sca])

out_df.to_csv(args["output"], sep="\t", header=False, index=False)
