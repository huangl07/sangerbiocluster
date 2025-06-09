# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据annosummry,将vcf 中exon区域的点提取出来，并进行
"""

import sys

def translate(vcf,gtf,output):
    exon_dict = {}
    with open(gtf,"r")as g:
        i = 0
        for line in g:
            if i == 0:
                i+=1
                pass
            else:
                lines =line.strip().split("\t")
                if lines[2] == "exon":
                    chr = lines[0]
                    start = lines[3]
                    end = lines[4]
                    #print(start,end)
                    if chr not in exon_dict:
                        exon_dict[chr] = []
                        exon_dict[chr].append(start+":"+end)
                    else:
                        exon_dict[chr].append(start+":"+end)
    with open(output,"w") as o:
        with open(vcf,"r")as v:
            ind = 0
            for line in v:
                if "#" in line:
                    o.write(line)
                else:
                    exon_sign = False
                    chr = line.split("\t")[0]
                    pos = line.split("\t")[1]
                    if chr in exon_dict:
                        for value in exon_dict[chr]:
                            start = int(value.split(":")[0])
                            end = int(value.split(":")[1])
                            if start <= int(pos) <= end:
                                exon_sign = True
                            else:
                                pass
                        if exon_sign == True:
                            o.write(line)
                    else:
                        pass
                    #else:
                        #print(chr+":"+str(pos)+" versus "+ str(start)+":"+str(pos))
                    
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly two arguments: vcf, database,output and name. Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
