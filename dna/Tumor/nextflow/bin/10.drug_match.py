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
                gene_all = lines.strip().split("\t")[0]
                gene = gene_all.split(" ")[0]
                line = lines.strip().split("\t")
                out_line = []
                out_line.append(line[2])
                out_line.append(line[5])
                out_line.append(line[7])
                out_line.append(line[8])
                out_line.append(line[9])
                out_line.append(line[10])
                if gene not in gene_dict:
                    gene_dict[gene] = []
                    gene_dict[gene].append(out_line)
                else:
                    gene_dict[gene].append(out_line)
    with open(output,"w") as o:
        i = 0
        o.write("Symbol\tChrom\tPos\tRef\tAlt\tProtein_change\tDisease	Drug	Drug_Type	Drug_Directions	Evidence_Level	Significance\n")
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
                    classification = fields[7].split("|")[11]
                    if gene_id in gene_dict:
                        out_str = ""
                        for z in range(len(gene_dict[gene_id])):
                            out_str +=  gene_id + "\t" + "\t".join(line.strip().split("\t")[0:2]) + "\t" + "\t".join(line.strip().split("\t")[3:5])+ "\t"+ classification + "\t" + "\t".join(gene_dict[gene_id][z]) + "\n"
                        o.write(out_str)
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly two arguments: vcf, database and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
