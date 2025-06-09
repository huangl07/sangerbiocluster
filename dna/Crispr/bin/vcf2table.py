# !usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__: yuan.xu
# first_modified: 20220831
# last_modified: 20221208

"""
把vcf转成成pop.table格式
主要使用vcfpy包
"""
import argparse
import vcfpy
from operator import add
from functools import reduce
from typing import List


def table_header_extend(sample_list: List[str]) -> List[str]:
    """
    功能：把sample_list的表头扩充，例如：[B1, B2] => [B1.GT, B1.AD, B1.DP, B2.GT, B2.AD, B2.DP]
        ：param sample_list 一个需要修改的的List
        : return List[str]
    """
    new_sample_list = []
    for sample_name in sample_list:
        sample_name_GT = sample_name + '.GT'
        sample_name_AD = sample_name + '.AD'
        sample_name_DP = sample_name + '.DP'
        new_sample_list.append(sample_name_GT)
        new_sample_list.append(sample_name_AD)
        new_sample_list.append(sample_name_DP)
    return new_sample_list


def call_trans(record_calls: List['Call']) -> List[str]:
    """
    功能：把record.calls对象转化成列表
        : param call_trans Call对象列表
        : return List[str]
    """
    gt_list = [str(call.gt_bases[0]) + "/" + str(call.gt_bases[1]) for call in record_calls]
    ad_list = [",".join(list(map(str, call.data.get('AD')))) for call in record_calls]  # 以,分隔符合并AD
    dp_list = [str(call.data.get('DP')) for call in record_calls]
    call_detail_info = reduce(add, [list(i) for i in zip(gt_list, ad_list, dp_list)])
    return call_detail_info

def abstract_info(info_list: List[str]) -> str:
    """
    功能：在有多个ann注释的情况下，只取第一个ann结果, vcf snpeff注释的info字段信息提取
        : info_list info字符串列表
        : return str
    """
    one_anno_info = info_list[0]
    one_anno_info_list = [i if i != "" else "--" for i in one_anno_info.split("|")]
    new_ann_info = "|".join([one_anno_info_list[0], one_anno_info_list[1], one_anno_info_list[2], one_anno_info_list[3], one_anno_info_list[4],
                         one_anno_info_list[10]])
    return new_ann_info

def abstract_info_multi(info_list: List[str]) -> str:
    """
    功能：在有多个ann注释的情况下，取多个注释, vcf snpeff注释的info字段信息提取
        : info_list info字符串列表
        : return str
    """
    new_info_list = []
    for info in info_list:
        one_anno_info_list = [i if i != "" else "--" for i in info.split("|")]
        new_one_ann_info = "|".join(one_anno_info_list)
        new_info_list.append(new_one_ann_info)
    return ";".join(new_info_list)


def trans_type(input_str: str) -> str:
    """
    功能：把SNV字符串转化成SNP;INS和DEL字符串转化成INDEL
        : input_str 要转化的字符串
        : return str 转化成的字符串
    """
    if input_str == "SNV":
        output_str = "SNP"
    elif input_str == "INS" or input_str == "DEL" or input_str == "INDEL":
        output_str = "INDEL"
    else:
        raise Exception("字符串不是SNV,INS,DEL,INDEL中的一种")
    return output_str

def vcf_to_table(infile: str, outfile: str) -> None:
    """
    主函数
    """
    # 读取文件(支持.gz形式)
    reader = vcfpy.Reader.from_path(infile)
    # 表头
    with open(outfile, "w", encoding='utf-8') as f:
        header = ['CHROM', 'POS', 'Ref', 'Alt', "Vtype"] + table_header_extend(reader.header.samples.names) + ["ANN"]
        f.write("\t".join(header))
        f.write("\n")
        # 读取每一个Record，获得所需要的信息
        for record in reader:
            if record.FILTER[0] == 'PASS' and len(record.ALT) == 1:
                try:
                    line = [str(record.CHROM), str(record.POS), record.REF]
                    line += [",".join([alt.value for alt in record.ALT])]
                    line += [",".join([alt.type for alt in record.ALT])]
                    line += call_trans(record.calls)
                    line += [abstract_info_multi(record.INFO['ANN'])]
                    f.write("\t".join(line))
                    f.write("\n")
                except Exception as err:
                    print(record.ALT)
                    print("An exception happened: " + str(err))
            else:
                pass


parser = argparse.ArgumentParser(description="vcf格式转化成pop.table")
parser.add_argument("-i", "--infile", required=True)
parser.add_argument("-o", "--outfile", required=True)
args = vars(parser.parse_args())
vcf_to_table(args["infile"], args["outfile"])
