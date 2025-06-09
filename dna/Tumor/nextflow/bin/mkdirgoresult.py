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

#fq_list=[]
for line in lines[1:]:
    tmp = line.strip().split("\t")
    chrname = tmp[0]
    goid = tmp[7]
    if len (goid) !=0 and goid != "--":
        x = goid.replace(",", "; ")
        print(chrname+'\t'+x)
    else:
        print(chrname)
