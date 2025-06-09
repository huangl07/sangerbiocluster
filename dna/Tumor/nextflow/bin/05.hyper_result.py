# -*- coding:utf-8 -*-
# @Last-edit Time 2023/05/08
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:根据report_dir路径下所有report结果中msi结果提取出来，并且合并
"""

import sys
import os
def rank(value,software):
    if software =="msi":
        if float(value) < 10:
            rank = "MSS"
        elif float(value) > 30:
            rank = "MSI-H"
        else:
            rank = "MSI-L"
    elif software == "msi_pro":
        if float(value) < 10:
            rank = "MSS"
        elif float(value) > 30:
            rank = "MSI-H"
        else:
            rank = "MSI-L"
    elif software == "msi2":
        if float(value) < 20:
            rank = "MSS"
        else:
            rank = "MSI-H"
    else:
        exit("ERROR:Wrong Software: "+ sys.argv[3] +". Exiting...")
    return rank
def translate(report_dir,output_dir,software):
    path = report_dir
    file_list = [f for f in os.listdir(path) if f.endswith('.prefix')]
    with open(os.path.join(output_dir,"msi_result_all.xls"),"w")as o:
        with open(os.path.join(output_dir,"msi_result_all.txt"),"w")as t:
            o.write("Tumor\tControl\tTotal_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tMSI_state\n")
            t.write("Tumor\tControl\tTotal_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tMSI_state\n")
            for file in file_list:
                with open(os.path.join(path,file),"r")as f:
                    sample = file.split(".")[0]
                    tumor = sample.split("_")[0]
                    control = sample.split("_")[1]
                    for lines in f:
                        if lines.startswith("Total_Number_of_Sites"):
                            pass
                        else:
                            output_str = lines.strip().split("\t")
                            msi_state_value = output_str[2]
                            msi_rank = rank(msi_state_value,software)
                            o.write(str(tumor)+"\t"+ str(control) + "\t"+"\t".join(output_str)+"\t"+ msi_rank +"\n")
                            t.write(str(tumor)+"\t"+ str(control) + "\t"+"\t".join(output_str)+"\t"+ msi_rank +"\n")
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly three arguments: report_dir,output_dir and software . Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[3]+" result")
    if sys.argv[3] == "msi" or sys.argv[3] == "msi_pro" or sys.argv[3] == "msi2":
        translate(sys.argv[1], sys.argv[2],sys.argv[3])
    else:
        exit("ERROR:Wrong Software: "+ sys.argv[3] +". Exiting...")
if __name__ == "__main__":
    main()
