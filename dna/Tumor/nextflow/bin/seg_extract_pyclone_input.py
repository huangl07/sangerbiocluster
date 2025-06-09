# Imports necessary modules
import csv, os

# Iterates through the files in the directory
with open('D:/236/cnl/sequenzaR包结果/Test_segments.txt') as file:
    sample_data={}
    sample_name = os.path.basename('D:/236/cnl/sequenzaR包结果/Test_segments.txt')
    sample_name = sample_name.split("_")[0]
    i = 0
    for lines in file:
        if i == 0:
            i = i+1
            pass
        else:
            line = lines.strip().split("\t")
            chr=line[0].split("chr")
            line_data = [chr[1], line[1], line[2], line[11], line[10]]
            index = sample_name+":"+line[0]+":"+line[1]+":"+line[2]
            sample_data[index]=line_data
            if len(chr) == 1:
                pass
            elif int(line[10])==0:
                print("major>0")
                pass
            else:
                line_data = [chr[1], line[1], line[2], line[11],line[10]]
                index = sample_name+":"+line[0]+":"+line[1]+":"+line[2]
                sample_data[index]=line_data
# Writes data to CSV
with open("D:/236/cnl/"+sample_name+"_combine_1.txt",'w') as outfile:
    title = ['chr', 'start', 'end', 'minor_cn', 'major_cn']
    outfile.write("\t".join(title)+"\n") # Write header row first
    for index, data in sample_data.items():
        outfile.write("\t".join(data)+"\n")

with open("D:/236/cnl/Test_snp.txt") as combine:
    sample_data={}
    sample_name = os.path.basename('D:/236/cnl/Test_snp.txt')
    sample_name = sample_name.split("_")[0]
    i = 0
    rsid=1
    for lines in combine:
        if i == 0:
            i = i+1
            pass
        else:
            line = lines.strip().split("\t")
            chr=line[0].split("chr")
            if len(chr) == 1:
                pass
            else:
                line_data = [chr[1], line[1], line[1], str("rs")+str(rsid),str(int(line[4])+int(line[5])),str(int(line[8])+int(line[9]))]
                index = sample_name+":"+line[0]+":"+line[1]+":"+line[1]
                sample_data[index]=line_data
                rsid+=1
# Writes data to CSV
with open("D:/236/cnl/"+sample_name+"_combine_2.txt",'w') as outfile2:
    title = ['chr','start', 'end',"rs", 'ref_counts', 'var_counts']
    outfile2.write("\t".join(title)+"\n") # Write header row first
    for index, data in sample_data.items():
        outfile2.write("\t".join(data)+"\n")
