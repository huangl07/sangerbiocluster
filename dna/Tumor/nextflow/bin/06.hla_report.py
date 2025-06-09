# -*- coding:utf-8 -*-
# @Last-edit Time 2023/05/05
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据report_dir路径下所有report结果中新抗原结果提取出来，并且合并
"""

import sys
import os

def translate(report_dir,output_dir):
    path = report_dir
    file_list = [f for f in os.listdir(path) if f.endswith('.report')]
    with open(os.path.join(output_dir,"MHC_result_all.xls"),"w")as o:
        with open(os.path.join(output_dir,"MHC_result_all.txt"),"w")as t:
            o.write("Sample\tNeoantigen\tType\tFeature\tExon\n")
            t.write("Sample\tNeoantigen\tType\tFeature\tExon\n")
            for file in file_list:
                with open(os.path.join(path,file),"r")as f:
                    sample = file.split(".")[0]
                    for lines in f:
                        if lines.startswith("HLA gene : "):
                            gene_name = lines.strip().split("HLA gene : ")[1]
                        elif lines.startswith("[Type 1]"):
                            output_str = lines.strip().split("\t")[1:]
                            o.write(str(sample)+"\t"+ gene_name + "\t"+"Type 1\t"+"\t".join(output_str)+"\n")
                            t.write(str(sample)+"\t"+ gene_name + "\t"+"Type 1\t"+"\t".join(output_str)+"\n")
                        
                        elif lines.startswith("[Type 2]"):
                            output_str = lines.strip().split("\t")[1:]
                            o.write(str(sample)+"\t"+ gene_name + "\t"+"Type 2\t"+"\t".join(output_str)+"\n")
                            t.write(str(sample)+"\t"+ gene_name + "\t"+"Type 2\t"+"\t".join(output_str)+"\n")
                        else:
                            pass
def main():
    if len(sys.argv) != 3:
        exit("ERROR: This program accepts exactly two arguments: vcf, database and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
