# -*- coding: utf-8 -*-
#!usr/bin/python3
from json.encoder import py_encode_basestring_ascii
import re
from sys import argv
import pandas as pd
import sys


file = open(sys.argv[1], encoding="utf-8", errors="ignore")
lines=file.readlines()
dict_CAPS_dCAPS={}
f = open(sys.argv[2], "a")
f.write("chr\tposition\tptype\tenzyme\tbases\tproduct_size\tTM_left\tGCcontent_left\tprimer_seq_left\tTM_right\tGCcontent_right\tprimer_seq_right\n")
for i in range(0,len(lines)):
    if(lines[i].startswith("index") and lines[i+1].strip() != ""):
        n=1
        t=0
        read1=lines[i+n]
        read2=lines[i+n+1]
        while read1.strip()!= "" and read2.strip() != "":
            enzyme=read1.split("-")[3].split(",")[0]
            ptype=read1.split("-")[2]
            pos=read1.split("-")[0]+"-"+read1.split("-")[1]
            chr = read1.split("-")[0]
            position =read1.split("-")[1]
            bases = read1.split("-")[4]
            product_size = read1.split("\t")[1]
            TM_left = read1.split("\t")[8]
            GCcontent_left = read1.split("\t")[9]
            primer_seq_left = read1.split("\t")[14]
            TM_right = read2.split("\t")[8]
            GCcontent_right = read2.split("\t")[9]
            primer_seq_right = read2.split("\t")[14]
            info = chr + "\t" + position + "\t" + ptype + "\t" + enzyme + "\t" + bases + "\t" + product_size + "\t" + TM_left + "\t" + GCcontent_left + "\t" + primer_seq_left + "\t" + TM_right + "\t" + GCcontent_right + "\t" + primer_seq_right

            if t == 0:
                dict_CAPS_dCAPS[pos]={}
                dict_CAPS_dCAPS[pos]["dCAPS"]={}
                dict_CAPS_dCAPS[pos]["CAPS"]={}
                t=1
            if ptype =="dCAPS" and len(dict_CAPS_dCAPS[pos]["dCAPS"].keys()) <3 and enzyme not in dict_CAPS_dCAPS[pos]["dCAPS"].keys():
                dict_CAPS_dCAPS[pos][ptype][enzyme] = info
            if ptype =="CAPS" and len(dict_CAPS_dCAPS[pos]["CAPS"].keys()) <3 and enzyme not in dict_CAPS_dCAPS[pos]["CAPS"].keys():
                dict_CAPS_dCAPS[pos][ptype][enzyme] = info
            i=i+2
            read1=lines[i+n]
            read2=lines[i+n+1]
        for v in dict_CAPS_dCAPS[pos]["dCAPS"].values():
            f.write(f"{v}\n")
        for v in dict_CAPS_dCAPS[pos]["CAPS"].values():
            f.write(f"{v}\n")
    else:
        pass
f.close()


file1 = open(sys.argv[3], encoding="utf-8", errors="ignore")
lines1=file1.readlines()
f1 = open(sys.argv[4], "a")
f1.write("chr\tposition\tproduct_size\ttype_1\tptype_1\tTM_1\tGCcontent_1\tprimer_seq_1\ttype_2\tptype_2\tTM_2\tGCcontent_2\tprimer_seq_2\ttype_3\tptype_3\tTM_3\tGCcontent_3\tprimer_seq_3\n")
dict_KASP={}
for i in range(0,len(lines1)):
    if(lines1[i].startswith("index") and lines1[i+1].strip()!= ""):
        read1=lines1[i+1]
        read2=lines1[i+2]
        read3=lines1[i+3]
        if read1.strip()!= "" and read2.strip() != "" and read3.strip() != "":
            pos = read1.split("-")[0] + "-" + read1.split("-")[1]
            chr = read1.split("-")[0]
            position = read1.split("-")[1]
            product_size = read1.split("\t")[1]
            type_1 = read1.split("\t")[2]
            ptype_1 = read1.split("-")[4]+"-"+read1.split("\t")[0].split("-")[5]
            TM_1 = read1.split("\t")[8]
            GCcontent_1= read1.split("\t")[9]
            primer_seq_1 = read1.split("\t")[14]
            type_2 =read2.split("\t")[2]
            ptype_2 =read2.split("-")[4]+"-"+read2.split("\t")[0].split("-")[5]
            TM_2 =read2.split("\t")[8]
            GCcontent_2 =read2.split("\t")[9]
            primer_seq_2 =read2.split("\t")[14]
            type_3 =read3.split("\t")[2]
            ptype_3 =read3.split("\t")[0].split("-")[4]
            TM_3 =read3.split("\t")[8]
            GCcontent_3 =read3.split("\t")[9]
            primer_seq_3 =read3.split("\t")[14]
            info = chr+"\t"+position+"\t"+product_size+"\t"+type_1+"\t"+ptype_1+"\t"+TM_1+"\t"+GCcontent_1+"\t"+primer_seq_1+"\t"+type_2+"\t"+ptype_2+"\t"+TM_2+"\t"+GCcontent_2+"\t"+primer_seq_2+"\t"+type_3+"\t"+ptype_3+"\t"+TM_3+"\t"+GCcontent_3+"\t"+primer_seq_3
            if not pos in dict_KASP.keys():
                dict_KASP.update({pos:info})
            i=i+3
for v in dict_KASP.values():
    f1.write(f"{v}\n")
f1.close()

#python3 capskasp_trimlist.py Potential_CAPS_primers.tsv 11 Potential_KASP_primers.tsv 111
