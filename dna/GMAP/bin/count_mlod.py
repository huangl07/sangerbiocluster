#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import argparse
from concurrent import futures
import pandas as pd
import time
def abstract_input(chunk):
    chr1 = chunk[0].str.split('-', n=1,expand=True)[0]
    chr2 = chunk[1].str.split('-', n=1,expand=True)[0]
    new_dataframe = pd.DataFrame({'Loci_1' : chr1, 'Loci_2' : chr2, 'MLOD' : round(chunk[2])})
    new_dataframe = new_dataframe.loc[chr1 != chr2]
    new_dataframe = exchange_Loci(new_dataframe)
    return new_dataframe

def all_chunk_append(iteror, thread):
    all_chunk_result = pd.DataFrame()
    with futures.ProcessPoolExecutor(thread) as pool:
        for chunk_result in pool.map(abstract_input, iteror):
            all_chunk_result = pd.concat([all_chunk_result,chunk_result])
    return all_chunk_result

###df Loci_1 Loci_2
def exchange_Loci(df):
    df['tmp'] = df['Loci_1'] <= df['Loci_2']
    df_f = df.loc[df['tmp'] == False]
    df_f = df_f.loc[:,['Loci_2','Loci_1','MLOD', 'tmp']]
    df_f.rename(columns={'Loci_2': 'Loci_1', 'Loci_1': 'Loci_2'}, inplace=True)
    df_t = df.loc[df['tmp'] == True]
    new_df = pd.concat([df_t, df_f])
    del new_df['tmp']
    return new_df

def max_MLOD_group(df, outfile1):
    summarize_df = df.groupby(["Loci_1", "Loci_2"])['MLOD'].max()
    chr_interval = df.groupby(["Loci_1", "Loci_2"])['MLOD'].max().index.to_list()[0]
    max_MLOD_group_result = chr_interval + (str(round(df.groupby(["Loci_1", "Loci_2"])['MLOD'].max(), 2)),)
    summarize_df.to_csv(outfile1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="mlod统计")
    parser.add_argument("-input", type = str, required = True, help = u"输入的文件*.mlod")
    parser.add_argument("-chunksize", type = int, required=False, default = 10000000, help = u"读取文件分块大小,默认10000000行一块")
    parser.add_argument("-thread", type=int, required=False, default = 12, help = u"进程池进程数,默认12")
    parser.add_argument("-output", type=str, required=False, default = 'summarize_MLOD_mean.csv', help = u"分组MLOD值平均值统计")
    args = parser.parse_args()
    input_file = args.input

    reader = pd.read_csv(input_file, sep = "\t", header = None, chunksize = args.chunksize, comment='#')
    start = time.time()
    all_df = all_chunk_append(reader, args.thread)
    max_MLOD_group(all_df, args.output)
    end = time.time()
    print("The time is %f s" % (end - start))
