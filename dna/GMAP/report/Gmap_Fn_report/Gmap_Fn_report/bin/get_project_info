#!/usr/bin/env python
#coding=utf-8

import os
import time
import argparse
parser=argparse.ArgumentParser(prog='get_project_info',usage='%(prog)s [opthions] [value]',description='This program is used to get project info!')

parser.add_argument('-C','--project_contract',help='the project contract',metavar='')
parser.add_argument('-M','--project_number',help='the project number',metavar='')
parser.add_argument('-K','--project_customer',help='the project customer',metavar='')
parser.add_argument('-g','--project_genome_chinese',help='the project genome_chinese',metavar='')
parser.add_argument('-G','--project_genome',help='the project genome',metavar='')
parser.add_argument('-N','--project_samplenum',help='the number of samles ',metavar='')
parser.add_argument('-O','--outfile',help='the output file',metavar='')
argv=vars(parser.parse_args())

if argv['project_contract'] == None:
    raise Exception('You should provide the project contract!')
else:
    project_contract=argv['project_contract'].strip()

if argv['project_number'] == None:
    raise Exception('You should provide the project number!')
else:
    project_number=argv['project_number'].strip()

if argv['project_customer'] == None:
    raise Exception('You should provide the project customer!')
else:
    project_customer=argv['project_customer'].strip()

if argv['project_genome_chinese'] == None:
    raise Exception('You should provide the project genome_chinese!')
else:
    project_genome_chinese=argv['project_genome_chinese'].strip()

if argv['project_genome'] == None:
    raise Exception('You should provide the project genome!')
else:
    project_genome=argv['project_genome'].strip()

if argv['project_samplenum'] == None:
    raise Exception('You should provide the number of samples!')
else:
    project_samplenum=argv['project_samplenum'].strip()

if argv['outfile'] == None:
    raise Exception('You should provide the output file!')
else:
    outfile=argv['outfile'].strip()

project_time=time.strftime('%Y-%m-%d',time.localtime(time.time()))

out=open(outfile,'w')
sup_out = open(outfile +'.txt','w')
out.write('项目信息\t项目内容\n')
out.write('合同编号\t%s\n' %(project_contract))
out.write('项目编号\t%s\n' %(project_number))
out.write('物种信息\t%s/%s\n' %(project_genome_chinese,project_genome))
sup_out.write('%s\t%s\t%s\n' %(project_genome_chinese,project_genome,project_samplenum))
out.write('客户姓名\t%s\n' %(project_customer))
out.write('报告时间\t%s\n' %(project_time))
out.close()
