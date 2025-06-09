# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20220908
# last_modified: 20220908

import argparse
from itertools import islice
import csv
from collections import namedtuple

def abstract_all_transcripts_from_anno_summary(infile: str, outfile: str) -> None:
    """
    主函数：从anno.summary文件中获取所有的transcript列表
    """
    with open(infile, "r", encoding='utf-8') as f, open(outfile, "w", encoding='utf-8') as m:
        header = ["chr", "transcript_start", "transcript_end", "transcript_id", "gene_start", "gene_end", "gene_id"]
        m.write("\t".join(header))
        m.write("\n")
        f_csv = csv.reader(f, delimiter="\t")
        headers = next(f_csv)
        anno_summary_info_nt = namedtuple('anno_summary_info_nt', headers)
        for row in f_csv:
            each_row = anno_summary_info_nt(*row)
            gene_id, transcript_id_info, transcript_chr, transcript_start, transcript_end = each_row.GeneID.split(":")
            transcript_id, *_ = transcript_id_info.split("|")
            new_line_list = [transcript_chr, transcript_start, transcript_end, transcript_id, transcript_start, transcript_end, gene_id]
            m.write("\t".join(new_line_list))
            m.write("\n")


parser = argparse.ArgumentParser(description="从gff文件中提取全部的transcript")
parser.add_argument("-i", "--infile", required=True)
parser.add_argument("-o", "--outfile", required=True)
args = vars(parser.parse_args())
abstract_all_transcripts_from_anno_summary(args["infile"], args["outfile"])

#测试用例：python3 bsa_abstract_all_transcripts_from_anno_summary.py --infile /mnt/ilustre/users/yuan.xu/data/BSA/代码测试/测试文件/anno.summary --outfile /mnt/ilustre/users/yuan.xu/data/BSA/代码测试/测试文件/genes_all.list