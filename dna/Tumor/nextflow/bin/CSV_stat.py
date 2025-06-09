import os
# Open the input file and read it line by line
with open('D:/236/cnl/PD4_snp.txt') as f:
    i = 0 
    dict_var = {}
    #Read each line of the file one by one
    for line in f:
        if i== 0:
            i = i + 1
            pass
        # Split the line on the comma delimiter to get a list of the data
        else:
            datalist = line.strip().split("\t")
            chr = datalist[0]
            pos = datalist[1]
            # Get the 5th and 6th columns
            col5_sum = int(datalist[4]) + int(datalist[5])
            # Get the 9th and 10th columns
            col9_sum = int(datalist[8]) + int(datalist [9])
            cnv_num= col9_sum - col5_sum
            # Print the sums
            print("Sum of the 5th and 6th columns:",col5_sum)
            print("Sum of the 9th and 10th columns:", col9_sum)
            print("cnv variation:", cnv_num)
            index = str(chr) + str(pos)
            dict_var[index] = cnv_num

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