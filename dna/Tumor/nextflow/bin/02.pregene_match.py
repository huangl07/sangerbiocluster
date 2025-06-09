# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据txt 结果将变异位点相关基因提取出来，并且和oncogene数据库进行比较
"""

import sys

def translate(txt,database,output):
    gene_dict = {}
    with open (database,"r") as d:
        ih = 0
        for lines in d:
            if ih == 0:
                ih+=1
                pass
            else:
                gene = lines.strip().split("\t")[0]
                if lines.strip().split("\t")[10] == "" or lines.strip().split("\t")[10] == "-":
                    cgc = lines.strip().split("\t")[9]
                else:
                    cgc = lines.strip().split("\t")[10]
                gene_dict[gene]=cgc
    with open(output,"w") as o:
        i = 0
        o.write("Chromosome	Pos	Ref	Alt	Rs_id	MetaSVM_pred	FATHMM_pred	LRT_pred	PROVEAN_pred	MutationTaster_pred	MutationAssessor_pred	Polyphen2_HDIV_pred	Polyphen2_HVAR_pred	SIFT_pred	Deleterious_count	Cosmic_id	Protein_change	AA_change	AF_eas	Effect	Impact	Gene	Geneid	CGC_Cancers\n")
        with open(txt, "r") as f:
            for line in f:
                if i == 0:
                    i += 1
                    pass
                else:
                    fields = line.strip().split("\t")
                    gene_id = fields[21]
                    if gene_id in gene_dict:
                        out_str = "\t".join(line.strip().split("\t")) + "\t"  + gene_dict[gene_id] + "\n"
                        o.write(out_str)
                    else:
                        o.write( "\t".join(line.strip().split("\t")) + "\t"+"-"+"\n")
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly two arguments: txt, database and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
