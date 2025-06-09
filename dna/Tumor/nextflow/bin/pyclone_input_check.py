with open('D:/236/cnl/TEST_pyclone_input.tsv','r') as infile, open('D:/236/cnl/TEST_pyclone_input_check.tsv','w') as outfile:
    header_read = infile.readline()
    outfile.write(header_read) # write header to output file
    for eachline in infile: 
        line_info = eachline.split('\t')
        if int(line_info[5]) <= 0:
            pass
        else:
            outfile.write(eachline)
