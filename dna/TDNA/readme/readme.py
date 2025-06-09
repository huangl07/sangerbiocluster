#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __time__ = 202303

import os,sys,random,argparse
#Bin = os.path.dirname(__file__)

######################################
parser = argparse.ArgumentParser(description="Creating the readme file of result.")
parser.add_argument("-f", "--file",type=str,default = None,help="The detailed information of file")
parser.add_argument("-t", "--typeinfo",type=str,default = None,help="The type of result. Eg:BSA,GWAS,Genetic")
parser.add_argument("-n", "--outname",type=str,default = 'WGS',help="The name of output")
parser.add_argument("-b", "--bin",type=str,default = '~/app/database/wgs_v4_bed/report_base/readme',help="The path of header.html and tail.html")
###############################################
args = parser.parse_args()
file_info,typeinfo,outname = args.file,args.typeinfo,args.outname
Bin = args.bin
######################################

if file_info == None and typeinfo == None:
    sys.exit('Please input infomation file of result.')
result_info = open(file_info, encoding='UTF8')
project_name = "WGS"

body_html = open('body.html','w', encoding='UTF8')

folder_rank = {}
for lines in result_info:
    
    line = lines.strip().split('\t')
    #print(line)
    if line[1] =="folder_1" :
        type_line = line[1]
        folder_rank['folder_rank_1'] = line[0]+'-'+str(random.random())
        file_rank = folder_rank['folder_rank_1']
        body_html.write('''<tr data-tt-id='%s'><td><span class='folder'>%s</a></span></td><td>%s</td><td>folder</td></tr> \n\n'''%(folder_rank['folder_rank_1'],line[0],line[2]))
        file_path = ('./'+line[0]).split()
        same_folder = 'folder_1'
    elif "folder" in line[1] and line[1] !="folder_1": 
        num = int(line[1].split('_')[1])
        parent_id = folder_rank['folder_rank_' +  str(num - 1)]
        random_id = line[0]+'-'+str(random.random())
        body_html.write('''<tr data-tt-id='%s' data-tt-parent-id='%s'><td><span class='folder'>%s</td><td>%s</td><td>folder</td></tr>\n\n'''%(random_id,parent_id,line[0],line[2]))
        folder_rank['folder_rank_' +  str(num)] = random_id
        file_rank = random_id 
        if int(line[1].split("_")[-1]) > int(same_folder.split("_")[-1]):
            file_path.append(line[0])
            same_folder = line[1]
        elif int(line[1].split("_")[-1]) == int(same_folder.split("_")[-1]):
            file_path.pop()
            file_path.append(line[0])
            same_folder = line[1]
        else:
            n = int(same_folder.split("_")[-1]) - int(line[1].split("_")[-1]) + 1
            while n > 0 :
                file_path.pop()
                n = n - 1
            file_path.append(line[0])
            same_folder = line[1]
    elif line[1] == "file":
        file_path_tmp = os.path.join("/".join(file_path),line[0])
        body_html.write('''<tr data-tt-id='%s' data-tt-parent-id='%s'><td><span class='file'><a\nhref='%s'>%s</a></span></td><td>%s</td><td>file</td></tr>\n\n'''%(line[0]+'-'+str(random.random()),file_rank,file_path_tmp,line[0],line[2]))
body_html.close()


os.system(''' cat %s/header.html body.html %s/tail.html > %s.ReadME.html && rm body.html '''%(Bin,Bin,outname))

os.system('''sed -i 's/项目/%s/g' %s.ReadME.html'''%(project_name,outname))

