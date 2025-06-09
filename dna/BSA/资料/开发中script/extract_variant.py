# -*- encoding: utf-8 -*-
'''
@File    :   extract_variant.py
@Time    :   2022/12/06 16:47:01
@Author  :   jingwang
@Version :   1.0
@Contact :   jing.wang1@majorbio.com
'''

# here put the import lib
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='help')
parser.add_argument('--table', '-t', help='pop.table')
parser.add_argument('--region', '-r', help='pop.region')
parser.add_argument('--outfile', '-o', help='文件前缀名，region.variant,region.eff')
args = parser.parse_args()

df = pd.DataFrame()
table= pd.read_csv(args.table,sep="\t")

with open(args.region,"r") as f:
    next(f)
    for line in f.readlines():
        chr = line.split("\t")[0]
        pos1 = int(line.split("\t")[1])
        pos2 = int(line.split("\t")[2])
        df1 = table[(table['CHROM'] == chr) & (table['POS']>=pos2) & (table['POS'] <=pos2)]
        df = pd.concat([df,df1])
eff = df[df['Ann'].str.contains('')]        
eff = df[(df['Ann'].str.contains('MODERATE')) | (df['Ann'].str.contains("HIGH"))]

df.to_csv(args.outfile + ".variant",index=False,sep='\t')
eff.to_csv(args.outfile + ".eff",index=False,sep='\t')