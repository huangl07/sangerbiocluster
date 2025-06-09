# -*- coding:utf-8 -*-
# @Last-edit Time 2023/03/13
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据tdna结果进行汇总整理
"""

import sys
import os
def translate(csv_path,output,insert_name):
    with open(output,"w") as o:
        o.write("Sample\tType\tRef_chrom\tStart\tEnd\n")
        path = csv_path
        csv_files = [f for f in os.listdir(path) if f.endswith('.csv')]
        for file in csv_files:
            sample_name = file.split(".csv")[0]
            with open(os.path.join(path, file), 'r') as f:
                i = 0
                for line in f:
                    if i == 0:
                        i += 1
                        pass
                    else:
                        cols = line.strip().split(',')
                        if cols[1] != insert_name:
                            if cols[2] != '':
                                o.write(sample_name+"\t"+cols[0]+"\t"+cols[1]+"\t"+cols[2]+"\t"+cols[3]+"\n")
                            if cols[11] != '':
                                o.write(sample_name+"\t"+cols[0]+"\t"+cols[1]+"\t"+cols[11]+"\t"+cols[12]+"\n")

def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly two arguments: vcf, database and output . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2],sys.argv[3])

if __name__ == "__main__":
    main()
