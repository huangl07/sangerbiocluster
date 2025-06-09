# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据bcftools query 结果将snp/indel前后50bp提取出来,做出fasta,方便后续MHC统计
"""

import pysam
import sys
import re

one_letter ={'Val':'V', 'Ile':'I', 'Leu':'L', 'Glu':'E', 'Gln':'Q','Asp':'D', 'Asn':'N', 'His':'H', 'Trp':'W', 'Phe':'F', 'Tyr':'Y','Arg':'R', 'Lys':'K', 'Ser':'S', 'Thr':'T', 'Met':'M', 'Ala':'A','Gly':'G', 'Pro':'P', 'Cys':'C'}

def get_seq(your_fasta_filename,chr_id, start, end): 
    print(chr_id)
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

def translate(vcf,fa,output):
    with open(vcf) as seq_in:
        
        with open(output,"w")as o:
            i = 0
            for lines in seq_in:
                if "#" in lines:
                    pass
                else:
                    line = lines.strip().split('\t')
                    i = i + 1
                    info = line[7]
                    protein_list = info.split("|")
                    protein = protein_list[7].split(".")[0]
                    if protein == "":
                        pass
                    else:
                        #print("protein: " + protein)
                        pos_nb = protein_list[11]
                        if pos_nb == "":
                            pass
                        else:
                            #print("pos_nb: " + pos_nb)
                            pos = re.findall("\d+",pos_nb)
                            if pos =="":
                                pass
                            else:
                                #print(pos)
                                pos = int(pos[0])
                                alt_ref = pos_nb[-3:]
                                if str(alt_ref) in one_letter.keys():
                                    alt = one_letter[pos_nb[-3:]]
                                    start = pos - 6
                                    if start <= 0:
                                        start = 0
                                        rep = pos+1
                                    else:
                                        rep = 7
                                    end = start +12
                                    seq = get_seq(fa,protein,start,end)
                                    #print(seq,"\t",start,end)
                                    if len(seq) == 0:
                                        pass
                                    else:
                                        seq_out = sub(seq,rep-2,alt)

                                        if len(seq_out) >= 12:
                                            #print("join:"+ seq_out)
                                                o.write(">"+line[0]+"_"+line[1]+"\n")
                                                o.write(str(seq_out)+"\n")
                                        else:
                                            pass
                                else: 
                                    #print(alt_ref) 
                                    if alt_ref == "dup":
                                        #print("dup")
                                        alt = pos_nb.split(".")[1][:3]
                                        if alt not in one_letter:
                                            print("alt:"+alt)
                                            print("pos_nb: " + pos_nb)
                                            print(alt_ref)
                                            print ("attention!~!!!!!!!!!")
                                        else:
                                            alt = one_letter[alt]
                                            #print(alt)
                                            start = pos - 6
                                            if start <= 0:
                                                start = 0
                                                rep = pos+1
                                            else:
                                                rep = 7
                                            end = start +12
                                            seq = get_seq(fa,protein,start,end)
                                            seq_o = list(seq)
                                            #print(seq_o)
                                            seq_o.insert(rep-2,alt)
                                            #print("insert: "+ seq_out)
                                            seq_out = "".join(seq_o)
                                            #print("join:"+ seq_out)
                                            if len(seq_out) >= 12:
                                                #print("join:"+ seq_out)
                                                o.write(">"+line[0]+"_"+line[1]+"\n")
                                                o.write(str(seq_out)+"\n")
                                            else:
                                                pass
                                    elif alt_ref == "del":
                                        alt = pos_nb.split(".")[1][:3]
                                        if alt not in one_letter:
                                            print("alt:"+alt)
                                            print("pos_nb: " + pos_nb)
                                            print(alt_ref)
                                            print ("attention!~!!!!!!!!!")
                                        else:
                                            alt = one_letter[alt]
                                            #print("delelte:" + alt)
                                            start = pos - 6
                                            if start <= 0:
                                                start = 0
                                                rep = pos+1
                                            else:
                                                rep = 7
                                            end = start +12
                                            seq = get_seq(fa,protein,start,end)
                                            #print(seq)
                                            seq_out = list(seq)
                                            seq_out.pop(rep-2)
                                            #print(seq_out)
                                            seq_out = "".join(seq_out)
                                            if len(seq_out) >= 12:
                                            #print("join:"+ seq_out)
                                                o.write(">"+line[0]+"_"+line[1]+"\n")
                                                o.write(str(seq_out)+"\n")
                                            else:
                                                pass
                                    else:
                                        print("pos_nb: " + pos_nb)
                                        print(alt_ref)
                                        print ("attention!~!!!!!!!!!")


def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly three arguments: the table1,table2 and table.stat.out . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2],sys.argv[3])

if __name__ == "__main__":
    main()
