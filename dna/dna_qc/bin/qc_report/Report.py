#!/usr/bin/env python
#coding=utf-8
'''
Created on 2021-08-26
author: yangwanyun
'''

import os
import sys
import argparse
import time
parser=argparse.ArgumentParser(prog='Report',usage='%(prog)s [opthions] [value]',description='This program is used to generate report files!')

parser.add_argument('-R','--result',help='the project result file',metavar='')
parser.add_argument('-P','--project',help='the project info ')
argv=vars(parser.parse_args())

#argparse

if argv['result'] == None:
    raise Exception('You should provide the result file!')
else:
    result=argv['result'].strip()

if argv['project'] == None:
    project_contract="WGS测序"
    project_number="WGS测序"
    project_customer="-"

else:
    project_info=argv['project'].strip()
    input=open(project_info,'r').readlines()[0].strip().split("\t")
    project_contract=input[0]
    project_number=input[1]
    project_customer=input[2]

name="WGS测序"+"report"
#result_dir=result+'/workflow_results'
result_dir=os.path.abspath(result)
report_dir=os.getcwd()+"/"+name
script = os.path.realpath(os.path.abspath(__file__))
pipline=os.path.dirname(script)
rmd_dir=pipline+'/rmd/'
rmd_list=[]
rmd_list.append(rmd_dir+'title.rmd')
rmd_list.append(rmd_dir+'workflow.rmd')
rmd_str=' '.join(rmd_list)

Step1Sh=open('./report.sh','w')
Step1Sh.write('if [ -d "%s" ]; then\n\n' %(report_dir))
Step1Sh.write('rm -rf %s\n\n' %(report_dir))
Step1Sh.write('fi\n\n')

Step1Sh.write('mkdir -p %s\n\n' %(report_dir))
Step1Sh.write('cp -r %s/src %s\n\n' %(pipline,report_dir))
Step1Sh.write('cp -r %s/css %s\n\n' %(pipline,report_dir))
Step1Sh.write('cat %s > %s/report.rmd\n\n' %(rmd_str,report_dir))
Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format html --outfile %s/report_html.rmd\n\n' %(pipline,report_dir,report_dir))
#Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format pdf --outfile %s/report_pdf.rmd\n\n' %(pipline,report_dir,report_dir))

Step1Sh.write('mkdir -p %s/file\n\n' %(report_dir))
Step1Sh.write('%s/bin/get_project_info --project_contract %s --project_number %s --project_customer %s --outfile %s/file/project_info.xls\n\n' %(pipline,project_contract,project_number,project_customer,report_dir))

Step1Sh.write('cp %s/01.CleanData/fig/*clean.base.png   %s/file/\n\n' %(result_dir,report_dir))
Step1Sh.write('cp %s/01.CleanData/fig/*clean.qual.png  %s/file/\n\n' %(result_dir,report_dir))
Step1Sh.write('cp %s/02.QCreport/02_QC/QC_stat.xls  %s/file/\n\n' %(result_dir,report_dir))
Step1Sh.write('cp %s/02.QCreport/qc-report.xls  %s/file/\n\n' %(result_dir,report_dir))

Step1Sh.write('Rscript  %s/bin/rmarkdown --rmd  %s/report_html.rmd --format html --outfile %s/report_raw.html\n\n' %(pipline,report_dir,report_dir))
#Step1Sh.write('Rscript  %s/bin/rmarkdown --rmd  %s/report_pdf.rmd --format pdf --outfile %s/report_for_pdf_raw.html \n\n' %(pipline,report_dir,report_dir))
Step1Sh.write('grep -v -E "Warning in stringi|integer vector"  %s/report_raw.html > %s/%s.html\n\n' %(report_dir,report_dir,name))
#Step1Sh.write('grep -v -E "Warning in stringi|integer vector" %s/report_for_pdf_raw.html > %s/%s_report_for_pdf.html\n\n' %(report_dir,report_dir,name))
Step1Sh.write('if [ -s %s/%s.html ];then\n\n' %(report_dir,name))
Step1Sh.write('    rm  %s/report.rmd\n\n' %(report_dir))
Step1Sh.write('    rm  %s/report_raw.html\n\n' %(report_dir))
#Step1Sh.write('    rm %s/report_for_pdf_raw.html\n\n' %(report_dir))
Step1Sh.write('    rm  %s/report_html.rmd\n\n'%(report_dir) )
#Step1Sh.write('    rm  %s/report_pdf.rmd\n\n'%(report_dir))
Step1Sh.write('    rm -rf  %s/css\n\n'%(report_dir))
Step1Sh.write('    rm -rf  %s/src\n\n'%(report_dir))
Step1Sh.write('else\n\n')
Step1Sh.write('exit 1\n\n')
Step1Sh.write('fi\n')
Step1Sh.close()
os.system('sh report.sh')
