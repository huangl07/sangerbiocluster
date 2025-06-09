# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据bcftools query 结果将snp/indel前后50bp提取出来,做出fasta,方便后续MHC统计
"""

import pysam
import sys

def get_seq(your_fasta_filename,chr_id, start, end): 
    fasta_open = pysam.Fastafile(your_fasta_filename) 
    sub_seq =fasta_open.fetch(chr_id, start, end)
    fasta_open.close()
    return sub_seq
def sub(string,p,c):#替换字符串string中指定位置p的字符为c
    new = []
    for s in string:
        new.append(s)
    new[p] = str(c)
    return ''.join(new)

def translate(si,fa,output):
    out = []
    with open(si) as seq_in:
        index = 0
        for lines in seq_in:
            if index == 0:
                index +=1
                pass
            else:
                line = lines.strip().split('\t')
                chr = line[0]
                pos = int(line[1])
                alt = line[3]
                start = pos - 50
                if start <= 0:
                    start = 0
                    rep = pos+1
                else:
                    rep = 51
                end = pos + 50
                seq = get_seq(fa,chr,start,end)
                seq_out = sub(seq,rep,alt)
                out.append(seq_out)
    with open(output,"w") as o:        
        for i in range(len(out)):
            id = i + 1
            out_id = "sca"+str(id)
            o.write(">"+out_id+"\n")
            o.write(out[i]+"\n")

def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly three arguments: the table1,table2 and table.stat.out . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2],sys.argv[3])

if __name__ == "__main__":
    main()
