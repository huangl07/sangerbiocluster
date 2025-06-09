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
args = parser.parse_args()
#
file_info,typeinfo,outname = args.file,args.typeinfo,args.outname
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
if(typeinfo =="BSA" or typeinfo =="BSR"):
	os.system('''cp readme_dir*/01.vcf2table_readme.txt 01.vcf2table/ && cp readme_dir*/02.ED-slid_readme.txt 02.ED-slid && cp readme_dir*/03.index-slid_readme.txt 03.index-slid/ && cp readme_dir*/04.Gprime_readme.txt	04.Gprime/ && cp readme_dir*/05.index-loess_readme.txt 05.index-loess && cp readme_dir*/06.enrich_readme.txt 06.enrich/ ''')
if(typeinfo =="BSA_new" or typeinfo =="BSR"):
    os.system('''cp readme_dir*/01.vcf2table_readme.txt 01.vcf2table/ && cp readme_dir*/02.ED_readme.txt 02.ED-slid && cp readme_dir*/03.index_readme.txt 03.index-slid/ && cp readme_dir*/04.Gprime_readme.txt	04.Gprime/ && cp readme_dir*/05.index-loess_readme.txt 05.index-loess && cp readme_dir*/06.enrich_readme.txt 06.enrich/ ''')
if(typeinfo == "BSA_mutmap"):
    os.system('''cp readme_dir*/01.vcf2table_readme.txt ./report/Result/data_release/01.vcf2table/ \
    && cp readme_dir*/02.index_readme.txt ./report/Result/data_release/02.index/ \
    && cp readme_dir*/03.loess_readme.txt ./report/Result/data_release/03.loess \
    && cp readme_dir*/04.enrich_readme.txt ./report/Result/data_release/04.enrich/ ''')
if(typeinfo == "BSA_ridit"):
    os.system('''cp readme_dir*/01.vcf2table_readme.txt ./report/Result/data_release/01.vcf2table/ \
    && cp readme_dir*/02.ridit_readme.txt ./report/Result/data_release/02.ridit/ \
    && cp readme_dir*/03.enrich_readme.txt ./report/Result/data_release/03.enrich/ ''')
if(typeinfo =="GPS"):
	os.system('''cp readme_dir/01.vcf_readme.txt 01.vcf/ && cp readme_dir/02.ridit_readme.txt 02.ridit/ && cp readme_dir/03.enrich_readme.txt 03.enrich/ ''')
elif(typeinfo =="GWAS"):
	os.system('''cp readme_dir/01.genotype_readme.txt 01.genotype/ && cp readme_dir/02.structure_readme.txt 02.structure/ && cp readme_dir/03.pca_readme.txt	03.pca/ && cp readme_dir/04.iqtree_readme.txt	04.iqtree/ && cp readme_dir/05.LD_readme.txt 05.LD/ && cp readme_dir/06.gwas_rMVP_readme.txt 06.gwas_rMVP && cp readme_dir/07.rMVP_region_anno_readme.txt 07.rMVP_region_anno/ && cp readme_dir/08.enrich_readme.txt 08.enrich''')
elif(typeinfo =="GWAS_yxy"):
        os.system('''cp readme_dir/01.GWAS.readme.txt 01.GWAS/ && cp readme_dir/02.GLM_enrich.readme.txt 02.GLM_enrich/ ''')
elif(typeinfo =="Genetic_notreemix"):
	project_name = "遗传进化"
	os.system('''cp readme_dir/01.vcf_filter_readme.txt 01.vcf_filter/ && cp readme_dir/02.tree_readme.txt 02.tree/ && cp readme_dir/03.structure_readme.txt 03.structure/ && cp readme_dir/04.pca_readme.txt 04.pca/ && cp readme_dir/05.genetic_diversity_readme.txt 05.genetic_diversity && cp readme_dir/06.psmc_readme.txt 06.psmc/ && cp readme_dir/07.dxy_readme.txt 07.dxy/ && cp readme_dir/08.LD_readme.txt 08.LD/ && cp readme_dir/09.sweep_readme.txt 09.sweep/''')
elif(typeinfo =="Genetic_treemix"):
	project_name = "遗传进化"
	os.system('''cp readme_dir/01.vcf_filter_readme.txt 01.vcf_filter/ && cp readme_dir/02.tree_readme.txt 02.tree/ && cp readme_dir/03.structure_readme.txt 03.structure/ && cp readme_dir/04.pca_readme.txt 04.pca/ && cp readme_dir/05.genetic_diversity_readme.txt 05.genetic_diversity && cp readme_dir/06.psmc_readme.txt 06.psmc/ && cp readme_dir/07.treemix_readme.txt 07.treemix/ && cp readme_dir/08.LD_readme.txt 08.LD/ && cp readme_dir/09.sweep_readme.txt 09.sweep/''')
elif(typeinfo =="Genetic_population"):
    project_name = "群体结构"
    os.system('''cp readme_dir/01.vcf_filter_readme.txt 01.vcf_filter/ && cp readme_dir/02.tree_readme.txt 02.tree/ && cp readme_dir/03.structure_readme.txt 03.structure/ && cp readme_dir/04.pca_readme.txt 04.pca/''')
elif(typeinfo =="gmap_F1"):
    project_name = "遗传图谱"
    os.system('''cp readme_dir/01.genotype_readme.txt 01.genotype/ && cp readme_dir/02.gmap_readme.txt 02.gmap && cp readme_dir/03.gmap_evaluate_readme.txt 03.gmap_evaluate/ && cp readme_dir/04.qtl_readme.txt 04.qtl && cp readme_dir/05.anno_readme.txt 05.anno && cp readme_dir/06.enrich_readme.txt 06.enrich/''')
elif(typeinfo =="gmap_Fn"):
    project_name = "遗传图谱"
    os.system('''cp readme_dir/01.genotype_readme.txt 01.genotype/ && cp readme_dir/02.gmap_readme.txt 02.gmap && cp readme_dir/03.gmap_evaluate_readme.txt 03.gmap_evaluate/ && cp readme_dir/04.qtl_readme.txt 04.qtl && cp readme_dir/05.anno_readme.txt 05.anno && cp readme_dir/06.enrich_readme.txt 06.enrich/''')
os.system('''sed -i 's/项目/%s/g' %s.ReadME.html'''%(project_name,outname))

