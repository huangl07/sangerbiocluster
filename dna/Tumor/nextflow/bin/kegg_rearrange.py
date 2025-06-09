#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np

def split_row(data, column):
  """拆分成行
  :param data: 原始数据
  :param column: 拆分的列名
  :type data: pandas.core.frame.DataFrame
  :type column: str
  """
  row_len = list(map(len, data[column].values))
  rows = []
  for i in data.columns:
    if i == column:
      row = np.concatenate(data[i].values)
    else:
      row = np.repeat(data[i].values, row_len)
    rows.append(row)
  return pd.DataFrame(np.dstack(tuple(rows))[0], columns=data.columns)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = r"将KEGG结果转化成列表")
    parser.add_argument("-filename", type = str, required = True, help = r"需要整理的KEGG结果")
    parser.add_argument("-outfile", type= str, required = True, help = r"结果输出文件")
    args = parser.parse_args()
    input_kegg_result = args.filename
    output_rearrange_result = args.outfile

    ###读入整理好的KEGG_result文件
    with open(input_kegg_result, "r") as f:
        lines = f.readlines()
        lines = [line.lstrip() for line in lines]

    ###把文件转换成以Query为key的字典
    KEGG_dict = {}
    for line in lines:
        if line.startswith("Query:"):
            id  = line.lstrip("Query:").strip()
            KEGG_dict[id] = ""
        else:
            KEGG_dict[id] = KEGG_dict[id] + line

    ###把字典用pandas转换成2列
    df = pd.DataFrame.from_dict(KEGG_dict, orient='index').reset_index()
    #改变两列名字分别为Query和Annotation
    df.rename(columns={'index': 'Query', 0: 'Annotation'}, inplace=True)
    #把Annotation列分割成KO和KEGG_pathway两列，目前有四列Query，Annotation，KO，KEGG_pathway
    df['KO'], df['KEGG_pathway'] = df['Annotation'].str.split('Pathway:').str[0], df['Annotation'].str.split('Pathway:').str[1]
    #KO列\t分割，取第一个K编号，KO：前缀
    df['KO'] = df['KO'].str.split('\t').str[1]
    #KEGG_pathway列根据\n分割成列表
    df.KEGG_pathway = df.KEGG_pathway.str.split('\n')
    #填充df缺少值
    df = df.fillna("-")
    #去除列表中的''元素
    def tmp(x):
        return [i for i in x if i]
    df.KEGG_pathway = df.KEGG_pathway.map(tmp)
    #拆分KEGG_pathway列成行
    df_arrange = split_row(df, "KEGG_pathway")
    #整理KEGG_pathway列并拆分成pathway_description和ko编号两列
    def tmp2(x):
        return x.strip()
    df_arrange.KEGG_pathway = df_arrange.KEGG_pathway.map(tmp2)
    df_arrange['pathway_description'], df_arrange['ko'] = df_arrange.KEGG_pathway.str.split('KEGG PATHWAY').str[0], df_arrange.KEGG_pathway.str.split('KEGG PATHWAY').str[1]
    df_arrange = df_arrange.fillna("-")
    df_arrange.pathway_description = df_arrange.pathway_description.map(tmp2)
    df_arrange.ko = df_arrange.ko.map(tmp2)

    ###输出整理完的结果
    df_arrange.to_csv(output_rearrange_result, columns=['Query', 'KO', 'pathway_description', 'ko'],index=False, sep='\t',quoting=3,escapechar=',')
