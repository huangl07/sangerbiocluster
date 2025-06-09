#!/usr/bin/env python
#coding=utf-8
'''
Created on 2022-05
author: jiawen.ma
'''

import os,glob
import sys
import argparse
import time
parser=argparse.ArgumentParser(prog='Report',usage='%(prog)s [opthions] [value]',description='This program is used to generate report files!')

parser.add_argument('-R','--result',help='the project result file',metavar='')
parser.add_argument('-E','--exclude',help='the exclude number',metavar='')
parser.add_argument('-P','--project',help='the project info ',metavar='')
argv=vars(parser.parse_args())

#argparse

if argv['result'] == None:
    raise Exception('You should provide the result file!')
else:
    result=argv['result'].strip()

if argv['exclude'] == None:
    exclude=str('0')
else:
    exclude=argv['exclude'].strip().split(',')

if argv['project'] == None:
    raise Exception('You should provide the project info file,Include project contract ,project number and customer name ,split by tap!')
else:
    project_info=argv['project'].strip()
    input=open(project_info,'r').readlines()[0].strip().split("\t")
    project_contract=input[0]
    project_number=input[1]
    project_customer=input[2]
    project_genome_chinese=input[3]
    project_genome=input[4]
    project_samplenum=input[5]

name=project_number+"_report"
#result_dir=result+'/workflow_results'
result_dir=os.path.abspath(result)
report_dir=os.getcwd()+"/"+name
script = os.path.realpath(os.path.abspath(__file__))
pipeline=os.path.dirname(script)
rmd_dir=pipeline+'/rmd/'
rmd_list=[]
rmd_list.append(rmd_dir+'title.rmd')
rmd_list.append(rmd_dir+'workflow.rmd')
for each in range(1,7):
    if str(each) not in exclude:
        rmd_list.append(rmd_dir+str(each)+'*.rmd')
rmd_list.append(rmd_dir+'appendix.rmd')
rmd_str=' '.join(rmd_list)


Step1Sh=open('./report.sh','w')
Step1Sh.write('if [ -d "%s" ]; then\n\n' %(report_dir))
Step1Sh.write('rm -rf %s\n\n' %(report_dir))
Step1Sh.write('fi\n\n')

Step1Sh.write('mkdir -p %s\n\n' %(report_dir))
Step1Sh.write('cp -r %s/src %s\n\n' %(pipeline,report_dir))
Step1Sh.write('cp -r %s/css %s\n\n' %(pipeline,report_dir))
Step1Sh.write('cat %s > %s/report.rmd\n\n' %(rmd_str,report_dir))
Step1Sh.write('''sed -i -e 's/genome_chinese/%s/g' -e 's/genome/%s/g' -e 's/sample_num/%s/g' %s/report.rmd \n\n'''%(project_genome_chinese,project_genome,project_samplenum,report_dir))
Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format html --outfile %s/report_html.rmd\n\n' %(pipeline,report_dir,report_dir))
#Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format pdf --outfile %s/report_pdf.rmd\n\n' %(pipeline,report_dir,report_dir))    

order=0
Step1Sh.write('mkdir -p %s/file\n\n' %(report_dir))
Step1Sh.write('%s/bin/get_project_info --project_contract %s --project_number %s --project_customer %s --project_genome_chinese %s --project_genome %s --project_samplenum %s --outfile %s/file/project_info.xls\n\n' %(pipeline,project_contract,project_number,project_customer,project_genome_chinese,project_genome,project_samplenum,report_dir))
filepath_01 = '%s/workflow_results/0%s.fastq_qc/cleandata_qc/atgc'%(result_dir,str(1))
sample_ID = os.path.split(glob.glob(os.path.join(filepath_01,'*xls'))[0])[1].split('.')[0]
Step1Sh.write('''sed -i 's/sample_Id/%s/g' %s/report_html.rmd\n\n''' %(sample_ID.split("-")[0],report_dir))
if '1' not in exclude:
    order=order+1
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/rawdata_qc/qc.xls   %s/file/rawdata.xls\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/cleandata_qc/atgc/*.atgc.xls  %s/file/\n\n'%(result_dir,str(order), report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/cleandata_qc/qual/*.qual.xls  %s/file/\n\n'%(result_dir,str(order), report_dir))
    Step1Sh.write('''sed -i '1d' %s/file/*.atgc.xls\n\n'''%(report_dir))
    Step1Sh.write('''cd %s/file && for j in `ls *atgc.xls`;do samples=`echo $j|sed 's/\.atgc\.xls//g'`;less $j |awk '{print $1"\t"$NF}'|sed '1d' > %s/file/$samples.qual.xls;Rscript %s/bin/ngsqc.r --base $samples.atgc.xls --qual $samples.qual.xls --key $samples --od %s/file/ ; done \n\n'''%(report_dir, report_dir,pipeline, report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/cleandata_qc/qc.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))

if '2' not in exclude:
    order=order+2
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/result.stat/Total.mapped.detail.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/insert/*.insert.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/coverage/*.coverage.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    #Step1Sh.write('cp %s/workflow_results/0%s.map_stat/depth/*.depth.xls  %s/file/\n\n' %(result_dir,str(order),sample_ID.split("-")[0],report_dir, sample_ID))
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/depth/*.depth.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('''cd %s/file && for j in `ls *insert.xls`;do sample=`echo $j|sed 's/\.insert\.xls//g'`; Rscript %s/bin/insert_size.R --i $j --o %s/file/$sample.insert.png ; done\n\n'''%(report_dir,pipeline,report_dir))
    #Step1Sh.write('Rscript %s/bin/insert_size.R --i %s/file/*.insert.xls --o %s/file/\n\n' %(pipeline,report_dir,report_dir))
    Step1Sh.write('less %s/file/Total.mapped.detail.xls |cut -f1-4 > %s/file/align_stat.xls \n\n' %(report_dir,report_dir))
    Step1Sh.write('''less %s/file/Total.mapped.detail.xls |awk -F "\t" '{print$1"\t"$8"\t"$9"\t"$6}' > %s/file/coverage_sample.xls \n\n'''%(report_dir,report_dir))
    #Step1Sh.write('''less %s/file/%s.coverage.xls |cut -f1|sort -u|grep -v 'sca'|sort -V > %s/file/chrlist \n\n''' %(report_dir,sample_ID,report_dir))
    Step1Sh.write('''less %s/file/%s.coverage.xls |cut -f1|sort -u|grep -v 'sca'|sort -V > %s/file/chrlist \n\n''' %(report_dir,sample_ID.split("-")[0],report_dir))
    Step1Sh.write('''cd %s/file && for j in `ls *depth.xls`;do samples=`echo $j|sed 's/\.depth\.xls//g'`;Rscript %s/bin/coverage_depth.R  --i $j --o %s/file/$samples.depth ; done\n\n''' %(report_dir,pipeline,report_dir))
    Step1Sh.write('''cd %s/file && for j in `ls *coverage.xls`;do samples=`echo $j|sed 's/\.coverage\.xls//g'`;Rscript %s/bin/genomeCoveragehorizontalArea.R --infile $j --idfile %s/file/chrlist --outfile %s/file/$samples.coverage --group.col 1 --x.col 2 --y.col 3 --x.lab Sequence-Position --y.lab AverageDepth-log2 --skip 0 --unit 100kb --log2 ; done\n\n''' %(report_dir,pipeline,report_dir,report_dir))


if '3' not in exclude:
    order=order+1
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/snp.stat  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/anno_stat/snp.stat  %s/file/snp_anno.xls \n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('less %s/file/snp.stat |cut -f1-7 > %s/file/snp_stat.xls \n\n' %(report_dir,report_dir))

if '4' not in exclude:
    order=4
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/indel.stat  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/anno_stat/indel.stat  %s/file/indel_anno.xls\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('less %s/file/indel.stat |cut -f1-7 > %s/file/indel_stat.xls \n\n' %(report_dir,report_dir))
    for i in range(0,int(project_samplenum)):
        sample_id_s = os.path.split(glob.glob(os.path.join(filepath_01,'*xls'))[i])[1].split('.')[0].split("-")[0]
        Step1Sh.write('''less %s/workflow_results/0%s.snp_indel/variant_stat/indel.len |grep %s > %s/file/%s.indel.len \n\n''' %(result_dir,str(order),sample_id_s,report_dir,sample_id_s))
        Step1Sh.write('Rscript %s/bin/indel_len.R  --i %s/file/%s.indel.len --o %s/file/\n\n' %(pipeline,report_dir,sample_id_s,report_dir))
    #Step1Sh.write('Rscript %s/bin/indel_len.R  --i %s/workflow_results/0%s.snp_indel/variant_stat/indel.len --o %s/workflow_results/0%s.snp_indel/variant_stat/\n\n' %(pipeline,result_dir,str(order),result_dir,str(order)))
    #Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/%s.indel.png  %s/file/\n\n' %(result_dir,str(order),sample_ID.split("-")[0],report_dir))

if '5' not in exclude:
    order=18
    Step1Sh.write('cp  %s/%s.filter/pop.table   %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('''cd %s/file &&less pop.table|cut -f1,5 > pop.stat && less chrlist |sed 's/chr//g' > chr_num.txt && cat chr_num.txt |while read line; do less pop.stat |awk '{if($1==k)print $0}' k=$line|cut -f2|sort|uniq -c |paste - - |awk '{print k"\t"$3"\t"$1}' k=$line >> chr_snp_indel_stat.xls; done && sed -i '1iChromosome ID\tSNP Number\tInDel Number' chr_snp_indel_stat.xls \n\n''' %(report_dir))
    Step1Sh.write('cp  %s/%s.filter/BSA/02.index-slid/pop.index.index.png   %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/BSA/03.Gprime/Gprime.index.png   %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/BSA/04.loess/loess.index.png   %s/file/\n\n' %(result_dir,str(order),report_dir))

if '6' not in exclude:
    order=18
    Step1Sh.write('cp  %s/pop.ED.region   %s/file/pop.ED.region.xls\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/ED/GO_result/ED_GOenrichment.png   %s/file/ED_go_enrich.png\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/ED/KEGG_result/ED_KEGGenrichment.png   %s/file/ED_kegg_enrich.png\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/Gprime/GO_result/Gprime_GOenrichment.png   %s/file/Gprime_go_enrich.png\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/Gprime/KEGG_result/Gprime_KEGGenrichment.png   %s/file/Gprime_kegg_enrich.png\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/index/GO_result/index_GOenrichment.png   %s/file/index_go_enrich.png\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp  %s/%s.filter/*/05.enrich/index/KEGG_result/index_KEGGenrichment.png   %s/file/index_kegg_enrich.png\n\n' %(result_dir,str(order),report_dir))
    

Step1Sh.write('cd %s/file \n\n'%(report_dir))
Step1Sh.write('ls *pdf |awk -F ".pdf" \'{print "convert "$0" "$1".png"}\'|sh \n\n')
Step1Sh.write('rm *pdf\n\n')

Step1Sh.write('Rscript  %s/bin/rmarkdown --rmd  %s/report_html.rmd --format html --outfile %s/report_raw.html\n\n' %(pipeline,report_dir,report_dir))
#Step1Sh.write('Rscript  %s/bin/rmarkdown --rmd  %s/report_pdf.rmd --format pdf --outfile %s/report_for_pdf_raw.html \n\n' %(pipeline,report_dir,report_dir))
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
