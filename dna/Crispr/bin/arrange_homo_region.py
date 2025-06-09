# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/08/02
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

"""
结果文件整理
"""

import argparse
import csv
from collections import namedtuple

def arrange_result(infile: str, outfile: str, bedfile: str) -> None:
    with open(infile, "r") as f, open(outfile, "w") as w, open(bedfile, "w") as m:
        f_csv = csv.reader(f, delimiter="\t")
        headers = next(f_csv)
        headers = [i.replace("#" , "").replace(" ", "_") for i in headers]
        result_nt = namedtuple('result_nt', headers)
        header_list = ["Number", "Bulge type", "crRNA", "DNA", "Chromosome", "Position", "Direction", "Mismatches", "Bulge_Size"]
        w.write("\t".join(header_list))
        w.write("\n")
        bed_line_list_bak = []
        for index, row in enumerate(f_csv, start=1):
            number = "No" + str(index)
            each_row = result_nt(*row)
            bulge_type = each_row.Bulge_type
            crrna = each_row.crRNA
            dna = each_row.DNA
            chromosome = each_row.Chromosome
            location = each_row.Position
            start_position = int(location) - 50
            if start_position < 1:
                start_position = 1
            start_position = str(start_position)
            end_position = str(int(location) + 50)
            direction = each_row.Direction
            mismatches = each_row.Mismatches
            bulge_size = each_row.Bulge_Size
            new_line_list = [number, bulge_type, crrna, dna, chromosome, location, direction, mismatches, bulge_size]
            bed_line_list = [chromosome, start_position, end_position]
            # if len(bed_line_list_bak) != 0 and bed_line_list_bak[:3] == bed_line_list[:3]:
            #     bed_line_list[0] = bed_line_list[0] + "_" + str(index)
            bed_line_list_bak = bed_line_list
            w.write("\t".join(new_line_list) + "\n")
            m.write("\t".join(bed_line_list) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="整理crispritz同源区域结果")
    parser.add_argument("-i", "--infile", help="crispritz同源区域结果", required=True)
    parser.add_argument("-o", "--outfile", help="整理完的结果文件", required=True)
    parser.add_argument("-b", "--bedfile", help="整理完的bed文件", required=True)
    args = vars(parser.parse_args())
    arrange_result(args["infile"], args["outfile"], args["bedfile"])
