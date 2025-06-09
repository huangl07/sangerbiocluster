file_input = open("D:/236/cnl/sequenzaR包结果/Test_segments.txt", "r")
file_output = open("D:/236/cnl/sequenzaR包结果/Test_segments2TCGA.txt", "w")

for line in file_input: 
    parts = line.strip().split("\t")
    print(parts)
    file_output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(parts[0], parts[1], parts[2], parts[3], parts[4]))
    
file_output.close()
