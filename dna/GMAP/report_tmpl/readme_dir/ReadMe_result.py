#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = jiawen.ma
# __time__ = 202206

import os,sys,random,argparse
Bin = os.path.dirname(__file__)

######################################
parser = argparse.ArgumentParser(description="Creating the readme file of result.")
parser.add_argument("-f", "--file",type=str,default = None,help="The detailed information of file")
parser.add_argument("-t", "--typeinfo",type=str,default = None,help="The type of result. Eg:BSA,GWAS,Genetic")
parser.add_argument("-n", "--outname",type=str,default = 'result',help="The name of output")
parser.add_argument("-r", "--readme_dir",type=str,help="The name of readme_dir")
parser.add_argument("-s", "--data_release",type=str,help="The name of data_release")
args = parser.parse_args()
#
file_info,typeinfo,outname,readme_dir, data_release = args.file,args.typeinfo,args.outname,args.readme_dir, args.data_release
######################################
if file_info == None and typeinfo == None:
    sys.exit('Please input infomation file of result.')
elif file_info != None:
    project_name = ''
    result_info = open(file_info, encoding='UTF8')
elif file_info == None and typeinfo != None:
    project_name = typeinfo
    outname = typeinfo
    filename = typeinfo + '_result.info.txt'
    result_info = open(os.path.join(Bin,'bin',filename), encoding='UTF8')

body_html = open('body.html','w', encoding='UTF8')

folder_rank = {}
for lines in result_info:
    line = lines.strip().split('\t')
    if line[1] =="folder_1" :
        type_line = line[1]
        folder_rank['folder_rank_1'] = line[0]+'-'+str(random.random())
        #folder_rank['folder_rank_1'] = line[0]
        file_rank = folder_rank['folder_rank_1']
        body_html.write('''<tr data-tt-id='%s'><td><span class='folder'>%s</a></span></td><td>%s</td><td>folder</td></tr> \n\n'''\
            %(folder_rank['folder_rank_1'],line[0],line[2]))
        file_path = ('./'+line[0]).split()
        same_folder = 'folder_1'
        # print(file_path)

    elif "folder" in line[1] and line[1] !="folder_1": 
        num = int(line[1].split('_')[1]) 
        parent_id = folder_rank['folder_rank_' +  str(num - 1)]
        random_id = line[0]+'-'+str(random.random())
        body_html.write('''<tr data-tt-id='%s' data-tt-parent-id='%s'><td><span class='folder'>%s</td><td>%s</td><td>folder</td></tr>\n\n'''\
            %(random_id,parent_id,line[0],line[2]))
        folder_rank['folder_rank_' +  str(num)] = random_id
        #folder_rank['folder_rank_' +  str(num)] = line[0]
        file_rank = random_id 
        # if(line[1] != same_folder):
        #     file_path.append(line[0])
        #     same_folder = line[1]
        # else:
        #     file_path.pop()
        #     file_path.append(line[0])
        #     same_folder = line[1]
        # print("/".join(file_path))
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
        #print(file_path_tmp )
        body_html.write('''<tr data-tt-id='%s' data-tt-parent-id='%s'><td><span class='file'><a\nhref='%s'>%s</a></span></td><td>%s</td><td>file</td></tr>\n\n'''\
            %(line[0]+'-'+str(random.random()),file_rank,file_path_tmp,line[0],line[2]))
body_html.close()
os.system(''' cat %s/bin/header.html body.html %s/bin/tail.html > %s.ReadME.html && rm body.html '''%(Bin,Bin,outname))
#os.system('''sed -i 's/项目/%s/g' %s.ReadME.html'''%(project_name,outname))
if(typeinfo =="GMAP"):
    project_name = "遗传图谱"
elif(typeinfo =="GMAP_qtl"):
    project_name = "遗传图谱"
os.system('''sed -i 's/项目/%s/g' %s.ReadME.html'''%(project_name,outname))

