# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20220907
# last_modified: 20220907

import argparse
import re
from itertools import islice
import csv
from collections import namedtuple

def abstract_anno(infile: str, outfile: str, go_anno: str, kegg_anno:str) -> None:
    """
    主函数：从anno.summary文件中获取GO和KEGG注释结果
    """
    with open(infile, "r", encoding='utf-8') as f, open(outfile, "w", encoding='utf-8') as m, \
        open(go_anno, "w", encoding='utf-8') as n, open(kegg_anno, "w", encoding='utf-8') as p:
            header = ["transcript_id", "ko_id", "ko_anno", "go_term", "go_anno"]
            kegg_header = ["Query", "pathway_description", "ko"]
            p.write("\t".join(kegg_header))
            p.write("\n")
            m.write("\t".join(header))
            m.write("\n")
            f_csv = csv.reader(f, delimiter="\t")
            headers = next(f_csv)
            anno_summary_info_nt = namedtuple('anno_summary_info_nt', headers)
            for row in f_csv:
                each_row = anno_summary_info_nt(*row)
                _, transcript_id_info, *_ = each_row.GeneID.split(":")
                transcript_id, *_ = transcript_id_info.split("|")
                new_line_list = [transcript_id, each_row.KoID, each_row.Koanno, each_row.GOTERM, each_row.GOANNO]
                m.write("\t".join(new_line_list))
                m.write("\n")
                new_go_term = each_row.GOTERM.replace(",", "; ")
                go_new_line_list = [transcript_id, new_go_term]
                n.write("\t".join(go_new_line_list))
                n.write("\n")
                ko_id_list = each_row.KoID.split(",")
                ko_anno_list = each_row.Koanno.split(":")
                if len(ko_id_list) >= 2:
                    for i in range(0, len(ko_id_list)):
                        kegg_new_line_list = [transcript_id, ko_anno_list[i], ko_id_list[i]]
                        p.write("\t".join(kegg_new_line_list))
                        p.write("\n")
                elif ko_id_list==[""]:
                    kegg_new_line_list = [transcript_id, "-", "-"]
                    p.write("\t".join(kegg_new_line_list))
                    p.write("\n")


parser = argparse.ArgumentParser(description="从anno.summary文件中获取GO和KEGG注释结果")
parser.add_argument("-i", "--infile", required=True)
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("-g", "--go_anno", required=True)
parser.add_argument("-k", "--kegg_anno", required=True)
args = vars(parser.parse_args())
abstract_anno(args["infile"], args["outfile"], args["go_anno"], args["kegg_anno"])

#测试用例：python3 bsa_abstract_anno.py --infile /mnt/ilustre/users/yuan.xu/data/BSA/代码测试/测试文件/anno.summary --outfile /mnt/ilustre/users/yuan.xu/data/BSA/代码测试/测试文件/go_kegg_info.list