# -*- coding: utf-8 -*-

import sys
import os
import gzip

#python t_get_seq.py input.vcf output.txt

VCF=sys.argv[1]
SED = sys.argv[2]

# def un_gz(file_name):
#     """ungz zip file"""
#     f_name = file_name.replace(".gz", "")
#     #获取文件的名称，去掉
#     g_file = gzip.GzipFile(file_name)
#     #创建gzip对象
#     open(f_name, "w+").write(g_file.read())
#     #gzip对象用read()打开后，写入open()建立的文件里。
#     g_file.close()
#     #关闭gzip对象

# unfile =un_gz(VCF)


# def read_gz_file(path):
#     if os.path.exists(path):
#         with gzip.open(path, 'r') as pf:
#             for line in pf:
#                 yield line
#     else:
#         print('the path [{}] is not exist!'.format(path))

# con = read_gz_file(VCF)
# if getattr(con, '__iter__', None):
with gzip.open(VCF, "rb") as f:

    with open(SED, "w") as m:
        for line in f:
            s = line.decode()
            if s.startswith('#'):
                pass;
            else:
                line_list = s.strip().split("\t")
                alle=line_list[4].split(",")
                if len(alle) > 1:
                    continue
                m.write("{}\t{}\t{}\t{}\n".format(line_list[0],line_list[1],line_list[3][-1],line_list[4][-1]))

sys.exit()


