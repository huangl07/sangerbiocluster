# -*- coding: utf-8 -*-
#!usr/bin/python3
#dandan.zhang

#把有参变异检测出来的anno.summary文件变成现在bsa流程中可以用的文件
#随便吧  一模一样算了

import re
from sys import argv
import sys

f1 = open(sys.argv[1], encoding="utf-8", errors="ignore")
lines = f1.readlines()

print("##ko"+ '\t'+"KEGG Orthology")
print("##Method: BLAST Options: evalue <= 1e-05; rank <= 5")
print("##Summary:      34211 succeed, 49190 fail" + '\n')
print("#Query  KO ID|KO name|Hyperlink")

#fq_list=[]
for line in lines:
    tmp = line.strip().split("\t")
    chrname = tmp[0]
    koid = tmp[-1]
    if "K" in koid:
        koid = koid
        koname = tmp[4]
        if 'GN=' in koname:
            tmp1 = koname.strip().split(" ")
            for i in range(0,len(tmp1)):
                if 'GN=' in tmp1[i]:
                    koannno = tmp1[i].split("GN=")[1]
                    print(chrname+'\t'+koid+"|"+koannno+"|http://www.genome.jp/dbget-bin/www_bget?ko:"+koid)
                else:
                    pass
    else:
        koid = "None"
        print(chrname+'\t'+koid)
print("\n"+"\n")

f1 = open(sys.argv[1], encoding="utf-8", errors="ignore")
lines = f1.readlines()
print("--------------------"+'\n'+'\n')
#fq_list=[]
for line in lines:
    tmp = line.strip().split("\t")
    chrname = tmp[0]
    KOID = tmp[-1]
    koid = tmp[5]
    if "," in koid:
        koid = koid.strip().split(",")[0]
    koanno = tmp[6]
    if ":" in koanno:
        koanno = koanno.strip().split(":")[0]
    koname = tmp[4]
    if "K" in KOID:
        if 'GN=' in koname:
            tmp1 = koname.strip().split(" ")
            for i in range(0,len(tmp1)):
                if 'GN=' in tmp1[i]:
                    konameRael = tmp1[i].split("GN=")[1]
                    #if koanno is not None:
                    if len (koanno) !=0:
                        print("////"+ '\n' + "Query:" + "\t"+ chrname +"\n"+ "KO:" +'\t'+KOID+'\t'+konameRael+ "\n"+"Pathway:"+"\t"+koanno+'\t'+"KEGG PATHWAY"+'\t'+koid)
                    else:
                        print("////"+ '\n' + "Query:" + "\t"+ chrname +"\n"+ "KO:" +'\t'+KOID+'\t'+konameRael)
                else:
                    pass
    if "--" in KOID:
        KOID = "None"
        print("////"+ '\n' + "Query:" + "\t"+ chrname)

print("////")

