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
parser.add_argument('-T','--filetype',default='True',help='The filetype is "True" if workflow_results is new version, otherwise filetype is "False"')
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
result_dir=result+'/workflow_results'
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
Step1Sh.write('''sed -i -e 's/genome_chinese/%s/g' -e 's/genome__latin/%s/g' -e 's/sample_num/%s/g' %s/report.rmd \n\n'''%(project_genome_chinese,project_genome,project_samplenum,report_dir))
Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format html --outfile %s/report_html.rmd\n\n' %(pipeline,report_dir,report_dir))
#Step1Sh.write('%s/bin/filter_rmd --rmd %s/report.rmd --format pdf --outfile %s/report_pdf.rmd\n\n' %(pipeline,report_dir,report_dir))    

#updata workflow/01.fastq_qc
if argv['filetype'] == 'True':
    os.system('''mkdir -p %s/workflow_results/01.fastq_qc/cleandata_qc/{atgc,qual} && mkdir -p %s/workflow_results/01.fastq_qc/rawdata_qc/{atgc,qual} && cd %s/workflow_results/01.fastq_qc/qc_stat  && for j in *json; do newsampleID=`echo $j|cut -d '.' -f1`;perl %s/bin/fastp.pl -i $j -o $newsampleID; less ${newsampleID}.clean.atgcn|awk -F "\t" '{print $1"\t"$2*100"\t"$3*100"\t"$4*100"\t"$5*100"\t"$6*100}'| sed "1i#pos\tA\tT\tG\tC\tN" > ../cleandata_qc/atgc/${newsampleID}.atgc.xls ;sed "1i#pos\tAver" ${newsampleID}.clean.qual > ../cleandata_qc/qual/${newsampleID}.qual.xls;done && cat *stat|grep -v '^#'|awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$6*100"\t"$5*100}'|sed "1isample ID\tRaw Reads\tRaw Bases(bp)\tRaw GC(%%)\tRaw Q30(%%)" > ../rawdata_qc/qc.xls && cat *stat|grep -v '^#'|awk -F "\t" '{print $1"\t"$8"\t"$9"\t"$12*100"\t"$11*100}'|sed "1isample ID\tClean Reads\tClean Bases(bp)\tClean GC(%%)\tClean Q30(%%)" > ../cleandata_qc/qc.xls'''%(result_dir,result_dir,result_dir,pipeline))


order=0
Step1Sh.write('mkdir -p %s/file\n\n' %(report_dir))
Step1Sh.write('%s/bin/get_project_info --project_contract %s --project_number %s --project_customer %s --project_genome_chinese %s --project_genome %s --project_samplenum %s --outfile %s/file/project_info.xls\n\n' %(pipeline,project_contract,project_number,project_customer,project_genome_chinese,project_genome,project_samplenum,report_dir))
filepath_01 = '%s/workflow_results/0%s.fastq_qc/cleandata_qc/atgc'%(result_dir,str(1))
sample_ID = os.path.split(glob.glob(os.path.join(filepath_01,'*xls'))[0])[1].split('.')[0]
Step1Sh.write('''sed -i 's/sample_Id/%s/g' %s/report_html.rmd\n\n''' %(sample_ID.split("-")[0],report_dir))
Step1Sh.write('cp %s/project.info.xls   %s/file\n\n' %(result_dir,report_dir))
Step1Sh.write('cp %s/result.info.xls   %s/file\n\n' %(result_dir,report_dir))
#Adding project info
Step1Sh.write('''cat %s/file/result.info.xls|while read line;do key=`echo "$line"|cut -f1`;value=`echo "$line"|cut -f2`;sed -i "s/${key}/${value}/g"  %s/report_html.rmd;done \n\n'''%(report_dir,report_dir) )

if '1' not in exclude:
    order=order+1
    Step1Sh.write('''sed 's/-1//g' %s/workflow_results/0%s.fastq_qc/rawdata_qc/qc.xls  >  %s/file/rawdata.xls\n\n''' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/cleandata_qc/atgc/*.atgc.xls  %s/file/\n\n'%(result_dir,str(order), report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.fastq_qc/cleandata_qc/qual/*.qual.xls  %s/file/\n\n'%(result_dir,str(order), report_dir))
    Step1Sh.write('''sed 's/-1//g' %s/workflow_results/0%s.fastq_qc/cleandata_qc/qc.xls >  %s/file/qc.xls \n\n''' %(result_dir,str(order),report_dir))
    Step1Sh.write('''sed -i  's/#sampleID/Sample ID/g' %s/file/qc.xls && sed -i 's/#sampleID/Sample ID/g' %s/file/rawdata.xls &&sed -i '1d' %s/file/*.atgc.xls\n\n'''%(report_dir,report_dir,report_dir))
    Step1Sh.write('''cd %s/file && for j in `ls *.qual.xls`;do samples=`echo $j|sed 's/\.qual\.xls//g'`;less $j |awk '{print $1"\t"$NF}'|sed '1d' > %s/file/$samples.qual.sort.xls;Rscript %s/bin/ngsqc.r --base $samples.atgc.xls --qual $samples.qual.sort.xls --key $samples --od %s/file/ ; done \n\n'''%(report_dir, report_dir,pipeline, report_dir))

if '2' not in exclude:
    order=order+2
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/result.stat/Total.mapped.detail.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.map_stat/coverage/*.coverage.xls  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('''less %s/file/Total.mapped.detail.xls |cut -f1-4|sed 1d |sed "1isample ID\tMapped Ratio(%%)\tProperly Mapped(%%)\tDuplicate Ratio(%%)" > %s/file/align_stat.xls \n\n''' %(report_dir,report_dir))
    Step1Sh.write('''less %s/file/Total.mapped.detail.xls |awk -F "\t" '{print$1"\t"$8"\t"$9"\t"$6}' |sed 1d |sed '1iSample ID\tCoverage 1X(%%)\tCoverage 5X(%%)\tAverage Depth'> %s/file/coverage_sample.xls \n\n'''%(report_dir,report_dir))
    Step1Sh.write('''less  %s/chr.list |cut -f1|sort -u|sort -V > %s/file/chrlist \n\n''' %(result_dir,report_dir))
    #Step1Sh.write('''less %s/file/%s.coverage.xls |cut -f1|sort -u|grep -v 'sca'|sort -V > %s/file/chrlist \n\n''' %(report_dir,sample_ID,report_dir))
    #
    # 
    # Step1Sh.write('''less %s/file/%s.coverage.xls |cut -f1|grep -v '^#'|sort -u|grep -v 'sca'|sort -V > %s/file/chrlist \n\n''' %(report_dir,sample_ID.split("-")[0],report_dir))
    # 这行修改  如果有chr chrlist则为chr  如果没有  则为sca前20 uniq|head -n 20
    #Step1Sh.write('''less %s/file/%s.coverage.xls |cut -f1|grep -wE "sca1|sca2|sca3|sca4|sca5|sca6|sca7|sca8"|grep -v '^#'|sort -u|sort -V > %s/file/chrlist \n\n''' %(report_dir,sample_ID.split("-")[0],report_dir))
    
    Step1Sh.write('''cd %s/file && for j in `ls *coverage.xls`;do samples=`echo $j|sed 's/\.coverage\.xls//g'`;Rscript %s/bin/genomeCoveragehorizontalArea.R --infile $j --idfile %s/file/chrlist --outfile %s/file/$samples.coverage --group.col 1 --x.col 2 --y.col 5 --x.lab Sequence-Position --y.lab AverageDepth-log2 --skip 0 --unit 100kb --log2 --title.lab \"genome coverage on $samples\"; done\n\n''' %(report_dir,pipeline,report_dir,report_dir))
    Step1Sh.write('''sed -i 's/^#//g' %s/file/align_stat.xls  %s/file/coverage_sample.xls\n\n''' %(report_dir,report_dir))


if '3' not in exclude:
    order=order+1
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/snp.stat  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/anno_stat/snp.stat  %s/file/snp_anno.xls \n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('''awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' %s/file/snp_anno.xls |grep -v "MODERATE" |grep -v "HIGH"|grep -v "MODIFIER" |grep -v "LOW" |awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' > %s/file/snp_region_anno.xls \n\n''' %(report_dir,report_dir))
    Step1Sh.write('''awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' %s/file/snp_anno.xls |grep -E "sampleID|MODERATE|HIGH|MODIFIER|LOW" |awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' > %s/file/snp_effect_anno.xls \n\n''' %(report_dir,report_dir))

    Step1Sh.write('''less %s/file/snp.stat |cut -f1-7 |sed 's/SNPnumber/SNP Number/g'> %s/file/snp_stat.xls \n\n''' %(report_dir,report_dir))
    Step1Sh.write('''sed -i 's/^#//g' %s/file/snp_anno.xls \n\n''' %(report_dir))

if '4' not in exclude:
    order=4
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/indel.stat  %s/file/\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/anno_stat/indel.stat  %s/file/indel_anno.xls\n\n' %(result_dir,str(order),report_dir))
    Step1Sh.write('''awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' %s/file/indel_anno.xls |grep -v "MODERATE" |grep -v "HIGH"|grep -v "MODIFIER" |grep -v "LOW" |awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' > %s/file/indel_region_anno.xls \n\n''' %(report_dir,report_dir))
    Step1Sh.write('''awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' %s/file/indel_anno.xls |grep -E "sampleID|MODERATE|HIGH|MODIFIER|LOW" |awk '{ for (i=1; i<=NF; i++){if(NR==1){arr[i]=$i;}else{arr[i]=arr[i] "\t" $i}}}END{for (i=1; i<=NF; i++){print arr[i]}}' > %s/file/indel_effect_anno.xls \n\n''' %(report_dir,report_dir))
    Step1Sh.write('''less %s/file/indel.stat |cut -f1-5 |sed 's/Insert/Insertion/g' |sed 's/Delete/Deletion/g'> %s/file/indel_stat.xls \n\n''' %(report_dir,report_dir))
    #for i in range(0,int(project_samplenum)):
     #   sample_id_s = os.path.split(glob.glob(os.path.join(filepath_01,'*xls'))[i])[1].split('.')[0].split("-")[0]
      #  Step1Sh.write('''less %s/workflow_results/0%s.snp_indel/variant_stat/indel.len |grep %s > %s/file/%s.indel.len \n\n''' %(result_dir,str(order),sample_id_s,report_dir,sample_id_s))
       # Step1Sh.write('Rscript %s/bin/indel_len.R  --i %s/file/%s.indel.len --o %s/file/\n\n' %(pipeline,report_dir,sample_id_s,report_dir))
    #Step1Sh.write('Rscript %s/bin/indel_len.R  --i %s/workflow_results/0%s.snp_indel/variant_stat/indel.len --o %s/workflow_results/0%s.snp_indel/variant_stat/\n\n' %(pipeline,result_dir,str(order),result_dir,str(order)))
    #Step1Sh.write('cp %s/workflow_results/0%s.snp_indel/variant_stat/%s.indel.png  %s/file/\n\n' %(result_dir,str(order),sample_ID.split("-")[0],report_dir))

if '5' not in exclude:
    #order=1
    Step1Sh.write('cp  %s/Result_bsa/01.vcf2table/pop.table   %s/file/\n\n' %(result_dir,report_dir))
#    if argv['chrlist'] == None:
    Step1Sh.write('''cd %s/file &&less pop.table|cut -f1,5 > pop.stat && less chrlist > chr_num.txt && cat chr_num.txt |while read line; do SNP=$(less pop.stat |awk '{if($1==k)print $0}' k=$line|grep "SNP" |wc -l ) ; indel=$(less pop.stat |awk '{if($1==k)print $0}' k=$line|grep "INDEL" |wc -l );echo "$line $SNP $indel"  >> chr_snp_indel_stat.xls ;done &&  less -S chr_snp_indel_stat.xls |tr -s ' ' '\t' >chr_snp_indel_stat.xls1 && mv chr_snp_indel_stat.xls1 chr_snp_indel_stat.xls && sed -i '1iChromosome ID\tSNP Number\tInDel Number' chr_snp_indel_stat.xls\n\n''' %(report_dir))


    # Step1Sh.write('''cd %s/file &&less pop.table|cut -f1,5 > pop.stat && less chrlist > chr_num.txt && cat chr_num.txt |while read line; do less pop.stat |awk '{if($1==k)print $0}' k=$line|cut -f2|sort|uniq -c |paste - - |awk '{print k"\t"$3"\t"$1}' k=$line >> chr_snp_indel_stat.xls; done && sed -i '1iChromosome ID\tSNP Number\tInDel Number' chr_snp_indel_stat.xls \n\n''' %(report_dir))
#    else:
        # chrlist = argv['chrlist'].strip()
        # Step1Sh.write('''cd %s/file &&less pop.table|cut -f1,5 > pop.stat && less %s/../data/chr.list |cut -f1 > chr_num.txt && cat chr_num.txt |while read line; do less pop.stat |awk '{if($1==k)print $0}' k=$line|cut -f2|sort|uniq -c |paste - - |awk '{print k"\t"$3"\t"$1}' k=$line >> chr_snp_indel_stat.xls; done && sed -i '1iChromosome ID\tSNP Number\tInDel Number' chr_snp_indel_stat.xls \n\n''' %(report_dir,result_dir))

    #Step1Sh.write('''cd %s/file &&less pop.table|cut -f1,5 > pop.stat && less chrlist |sed 's/chr//g' > chr_num.txt && cat chr_num.txt |while read line; do less pop.stat |awk '{if($1==k)print $0}' k=$line|cut -f2|sort|uniq -c |paste - - |awk '{print k"\t"$3"\t"$1}' k=$line >> chr_snp_indel_stat.xls; done && sed -i '1iChromosome ID\tSNP Number\tInDel Number' chr_snp_indel_stat.xls \n\n''' %(report_dir))
    Step1Sh.write('cp  %s/Result_bsa/02.index-slid/pop.index.index.png   %s/file/1.index.index.png\n\n' %(result_dir,report_dir))
    # Step1Sh.write('cp  %s/Result_bsa/02.index-slid/pop.index.*.index.png   %s/file/\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/02.index-slid/pop.ED.index.png   %s/file/3.ED.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/04.loess/loess.index.png   %s/file/4.loess.index.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/03.Gprime/Gprime.index.png   %s/file/2.Gprime.png\n\n' %(result_dir,report_dir))
    # Step1Sh.write('cp  %s/Result_bsa/03.Gprime/Gprime.*.index.png   %s/file/\n\n' %(result_dir,report_dir))
    #Step1Sh.write('cp  %s/Result_bsa//04.loess/loess.index.png   %s/file/\n\n' %(result_dir,report_dir))

if '6' not in exclude:
    order=18
    Step1Sh.write('''if [ -s %s/Result_bsa/05.enrich/index/index.variant.stat ];then less -S %s/Result_bsa/05.enrich/index/index.variant.stat |awk '{print$1"_"$2"_"$3"\t""Index""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' |grep -v "Chromosome"> %s/file/index.variant.stat1;awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' %s/Result_bsa/02.index-slid/index.gene.total %s/file/index.variant.stat1 |awk '{if ($12 != "")print $7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' |sort -n -k 2 -k 3 |sed '1iMethond\tChrom\tPos Start\tPos End\t#Gene\t#GeneEff\tSNP\tInDel\t#Effsnp\t#EffInDel' > %s/file/1.index.pop.region.xls;fi \n\n''' %(result_dir,result_dir,report_dir,result_dir,report_dir,report_dir))

    Step1Sh.write('''if [ -s %s/Result_bsa/05.enrich/ED/ED.variant.stat ];then less -S %s/Result_bsa/05.enrich/ED/ED.variant.stat |awk '{print$1"_"$2"_"$3"\t""ED""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' |grep -v "Chromosome"> %s/file/ED.variant.stat1;awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' %s/Result_bsa/02.index-slid/ED.gene.total %s/file/ED.variant.stat1 |awk '{if ($12 != "")print $7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' |sort -n -k 2 -k 3 > %s/file/2.ED.pop.region.xls;fi \n\n''' %(result_dir,result_dir,report_dir,result_dir,report_dir,report_dir))

    Step1Sh.write('''if [ -s %s/Result_bsa/05.enrich/Gprime/Gprime.variant.stat ];then less -S %s/Result_bsa/05.enrich/Gprime/Gprime.variant.stat |awk '{print$1"_"$2"_"$3"\t""Gprime""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' |grep -v "Chromosome"> %s/file/Gprime.variant.stat1; awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' %s/Result_bsa/02.index-slid/Gprime.gene.total %s/file/Gprime.variant.stat1 |awk '{if ($12 != "")print $7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' |sort -n -k 2 -k 3> %s/file/3.Gprime.pop.region.xls ;fi\n\n''' %(result_dir,result_dir,report_dir,result_dir,report_dir,report_dir))

    Step1Sh.write('''if [ -s %s/Result_bsa/05.enrich/loess/loess.variant.stat ];then less -S %s/Result_bsa/05.enrich/loess/loess.variant.stat |awk '{print$1"_"$2"_"$3"\t""loess""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' |grep -v "Chromosome"> %s/file/loess.variant.stat1; awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' %s/Result_bsa/02.index-slid/loess.gene.total %s/file/loess.variant.stat1 |awk '{if ($12 != "")print $7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' |sort -n -k 2 -k 3> %s/file/4.loess.pop.region.xls;fi \n\n''' %(result_dir,result_dir,report_dir,result_dir,report_dir,report_dir))

    Step1Sh.write('''cat %s/file/*.pop.region.xls > %s/file/pop.region.xls\n\n''' %(report_dir,report_dir))

    Step1Sh.write('cp  %s/Result_bsa/05.enrich/index/GO_result/*.png   %s/file/1.index_go_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/index/KEGG_result/*.png   %s/file/1.index_kegg_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/ED/GO_result/*.png   %s/file/2.ED_go_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/ED/KEGG_result/*.png   %s/file/2.ED_kegg_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/Gprime/GO_result/*.png   %s/file/3.Gprime_go_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/Gprime/KEGG_result/*.png   %s/file/3.Gprime_kegg_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/loess/GO_result/*.png   %s/file/4.loess_go_enrich.png\n\n' %(result_dir,report_dir))
    Step1Sh.write('cp  %s/Result_bsa/05.enrich/loess/KEGG_result/*.png   %s/file/4.loess_kegg_enrich.png\n\n' %(result_dir,report_dir))

    # if os.path.exists(result_dir+ '/Result_bsa/02.index-slid/pop.index.region') and os.path.getsize(result_dir+ '/Result_bsa/02.index-slid/pop.index.region')!=0:
    #     Step1Sh.write("""less -S %s/Result_bsa/05.enrich/index/index.variant.stat |awk '{print$1"_"$2"_"$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > %s/file/index.variant.stat1 && awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$5"\t"$6"\t"$7"\t"$8}' %s/Result_bsa/02.index-slid/index.gene.total %s/file/index.variant.stat1 |cut -f 2-10 |sed 1d|sed '1 iChrom\tPos Start\tPos End\tGene\tGeneEff\tSNP\tInDel\tEffsnp\tEffInDel' > %s/file/pop.region.xls\n\n""" %(result_dir,report_dir,result_dir,report_dir,report_dir))
    #     #修改sca
    #     #Step1Sh.write('''less   %s/Result_bsa/02.index-slid/pop.index.region|grep -v -E 'KZ' > %s/file/pop.ED.region.xls\n\n''' %(result_dir,report_dir))

    #     Step1Sh.write('cp  %s/Result_bsa/05.enrich/index/GO_result/*.png   %s/file/index_go_enrich.png\n\n' %(result_dir,report_dir))
    #     Step1Sh.write('cp  %s/Result_bsa/05.enrich/index/KEGG_result/*.png   %s/file/index_kegg_enrich.png\n\n' %(result_dir,report_dir))

        # Step1Sh.write("""paste %s/Result_bsa/05.enrich/Gprime/Gprime.gene.stat %s/Result_bsa/05.enrich/Gprime/Gprime.variant.stat |cut -f 1,2,3,4,5,9,10,11,12 |awk '{print "Gprime""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'| grep -v "ChromosomeID"> %s/file/pop.Gprime.region.xls\n\n""" %(result_dir,result_dir,report_dir))
    
        # Step1Sh.write(""" cat %s/file/pop.index.region.xls %s/file/pop.Gprime.region.xls> %s/file/pop.region.xls\n\n""" %(report_dir,report_dir,report_dir))
    # else:
    #     Step1Sh.write("""paste %s/Result_bsa/05.enrich/Gprime/Gprime.gene.stat %s/Result_bsa/05.enrich/Gprime/Gprime.variant.stat |cut -f 1,2,3,4,5,9,10,11,12 |awk '{print "Gprime""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'| grep -v "ChromosomeID"> %s/file/pop.Gprime.region.xls && sed -i '1 i\Method\tChr\tpos_start\tpos_end\tGene\tGeneEff\t      SNP\tIndel\tEffsnp\tEffInDel' %s/file/pop.region.xls \n\n""" %(result_dir,result_dir,report_dir,report_dir))

    #     #Step1Sh.write('cp  %s/Result_bsa/05.enrich/ED/GO_result/ED_GOenrichment.png   %s/file/ED_go_enrich.png\n\n' %(result_dir,report_dir))
    #     #Step1Sh.write('cp  %s/Result_bsa/05.enrich/ED/KEGG_result/ED_KEGGenrichment.png   %s/file/ED_kegg_enrich.png\n\n' %(result_dir,report_dir))
    #     Step1Sh.write('cp  %s/Result_bsa/05.enrich/Gprime/GO_result/*.png   %s/file/Gprime_go_enrich.png\n\n' %(result_dir,report_dir))
    #     Step1Sh.write('cp  %s/Result_bsa/05.enrich/Gprime/KEGG_result/*.png   %s/file/Gprime_kegg_enrich.png\n\n' %(result_dir,report_dir))


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
if argv['filetype'] == 'True':
	Step1Sh.write('rm -r %s/workflow_results/01.fastq_qc/cleandata_qc && rm -r %s/workflow_results/01.fastq_qc/rawdata_qc \n\n'%(result_dir,result_dir))
Step1Sh.write('else\n\n')
Step1Sh.write('exit 1\n\n')
Step1Sh.write('fi\n')
Step1Sh.close()
os.system('sh report.sh')
