# -*- coding:utf-8 -*-
# @Last-edit Time 2023/04/04
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据gtf转换成exon_gene_list
"""
import sys
def translate(gtf,output):
    with open(output,"w") as o: 
        with open(gtf,"r")as g:         
            for lines in g:
                line = lines.strip().split("\t")
                if line[2]=="exon":
                    chr = line[0]
                    start = line[3]
                    end = line[4]
                    gene_name = line[8].split(";")[5]
                    gene = gene_name.split("\"")[1]
                    o.write(chr+"\t"+start+"\t"+end+"\t"+gene+"\n")
            else:
                pass
    
def main():
    if len(sys.argv) != 3:
        exit("ERROR: This program accepts exactly two arguments: gtf and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
