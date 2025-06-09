# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据vcf 结果将变异位点相关基因提取出来，并且和oncogene数据库进行比较
"""

import sys

def translate(vcf,database,output):
    gene_dict = {}
    with open (database,"r") as d:
        ih = 0
        for lines in d:
            if ih == 0:
                ih+=1
                pass
            else:
                gene = lines.strip().split("\t")[0]
                gene_dict[gene]=lines.strip().split("\t")[1:]
    with open(output,"w") as o:
        i = 0
        o.write("Symbol\tChrom\tPos\tRef\tAlt\tClassification\tIntOGen_Cancer_Type	IntOGen_Role	CGC_Tumour_Types	CGC_Cancer_Role	Com299_Pancaner_Frequency	Com435_Putative_DriveCategory	B125_Classification	B125_Core_pathway	B125_Process	SMG127_Pancan12_Freq\n")
        with open(vcf, "r") as f:
            for line in f:
                if i == 0:
                    i += 1
                    pass
                if line.startswith("#"):
                    # 如果是注释行，则直接跳过
                    pass
                else:
                    fields = line.strip().split("\t")
                    gene_id = fields[7].split("|")[4]
                    classification = fields[7].split("|")[2]
                    if gene_id in gene_dict:
                        out_str = gene_id + "\t" +"\t".join(line.strip().split("\t")[0:2]) + "\t" + "\t".join(line.strip().split("\t")[3:5]) + "\t" + classification + "\t" + "\t".join(gene_dict[gene_id]) + "\n"
                        o.write(out_str)
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly two arguments: vcf, database and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
