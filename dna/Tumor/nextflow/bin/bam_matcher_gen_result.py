# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/16
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据bam-matcher结果将Fraction of common和Conclusion结果提取出来
"""
import sys
import re
import os 
    
# Read the log file
def result(file):
    Conclusion = []
    Fraction = []
    Normal = []
    Tumor = []
    with open(file, 'r') as f: 
        # 逆序浏览文件
        lastLine = f.readlines()[-1].strip()
        # 输出最后一行的文本
        Conclusion.append(lastLine)

    with open(file, 'r') as f: 
        data = f.readlines() 
        
    # Extract the bam filenames and numbers 
        for line in data: 
            bam1 = re.findall(r"bam1\:.*\n", line) 
            if bam1:
                bam1_name = bam1[0].strip().split("/")[-1]
                bam1_name = bam1_name.split(".")[0]
                Tumor.append(bam1_name)
            bam2 = re.findall(r"bam2\:.*\n", line) 
            if bam2:
                bam2_name = bam2[0].strip().split("/")[-1]
                bam2_name = bam2_name.split(".")[0]
                Normal.append(bam2_name)
            fraction = re.findall(r"Fraction of common\:.*\n", line) 
            if fraction:
                out_fra = fraction[0].strip().split(":")[1]
                out_f = out_fra.split("(")[0]
                out_f =out_f.split(" ")[1]
                Fraction.append(out_f)
    re_str = str(Normal[0]) + "\t" + str(Tumor[0]) + "\t" + str(Fraction[0]) + "\t" + str(Conclusion[0]) + "\n"
    return re_str

def translate(bam_report_path,output):
    path = []
    for dirpath, _, filenames in os.walk(bam_report_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            # 筛选出.report后缀的文件，保存在列表中
            if os.path.splitext(filename)[1] == '.report':
                path.append(file_path)
    with open(output,"w")as o:
        #o.write("Normal\tTumor\tFraction\tConclusion\n")
        for pa in path:
            o.write(result(pa))

def main():
    if len(sys.argv) != 3:
        exit("ERROR: This program accepts exactly three arguments: the table1,table2 and table.stat.out . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
