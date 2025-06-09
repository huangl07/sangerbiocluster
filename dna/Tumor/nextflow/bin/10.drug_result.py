# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据vcf 结果将变异位点相关基因提取出来，并且和drug_target数据库进行比较
"""

import sys

def translate(vcf,database,output,name):
    gene_dict = {}
    variant_dict = {}
    with open (database,"r") as d:
        ih = 0
        for lines in d:
            if ih == 0:
                ih+=1
                pass
            else:
                variant_all =  lines.strip().split("\t")[9:12]
                variant =":".join(variant_all)
                line = lines.strip().split("\t")
                gene = line[1]
                out_line = line[1:]
                if gene not in gene_dict:
                    gene_dict[gene] = []
                    gene_dict[gene].append(out_line)
                else:
                    gene_dict[gene].append(out_line)
                if variant not in variant_dict:
                    variant_dict[variant] = []
                    variant_dict[variant].append(out_line)
                else:
                    variant_dict[variant].append(out_line)
    with open(output,"w") as o:
        i = 0
        o.write("Sample_Control\tChrom\tPos\tRef\tAlt\tMatching_Accuracy\tGene\tProtein_change\tDisease	Drug	Drug_Type	Drug_Directions	Evidence_Level	Significance\tDt_Chrom	Dt_Start	Dt_Stop	Dt_Ref	Dt_Alt	Evidence_Statement	Citation	Variant_origin	Representative_Transcript	Functional_Region\n")
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
                    variant_chr = fields[0]
                    variant_pos = fields[1]
                    gene_id = fields[7].split("|")[4] 
                    #classification = fields[7].split("|")[11]
                    for variant in variant_dict:
                        chr = variant.split(":")[0]
                        if variant_chr != chr:
                            pass
                        else:
                            start = variant.split(":")[1]
                            end = variant.split(":")[2]
                            if int(start)<=int(variant_pos) <= int(end):
                                out_str = ""
                                for q in range(len(variant_dict[variant])):
                                    out_str +=  name+"\t"+ "\t".join(line.strip().split("\t")[0:2]) + "\t" + "\t".join(line.strip().split("\t")[3:5])+ "\tVariants_Matched"+"\t" + "\t".join(variant_dict[variant][q]) + "\n"
                                o.write(out_str)
                            else:
                                if gene_id in gene_dict:
                                    out_str = ""
                                    for z in range(len(gene_dict[gene_id])):
                                        out_str +=  name+"\t"+ "\t".join(line.strip().split("\t")[0:2]) + "\t" + "\t".join(line.strip().split("\t")[3:5])+ "\tGenes_Matched"+ "\t" + "\t".join(gene_dict[gene_id][z]) + "\n"
                                    o.write(out_str)
def main():
    if len(sys.argv) != 5:
        exit("ERROR: This program accepts exactly two arguments: vcf, database,output and name. Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])

if __name__ == "__main__":
    main()
