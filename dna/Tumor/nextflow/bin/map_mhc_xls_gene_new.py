 # -*- coding:utf-8 -*-
# @Last-edit Time 2023/04/04
# @Author xiaoya.ye
# @mail xiaoya.ye@majorbio.com

"""
目的:map mhc_xls with gene
"""
import sys
def translate(gene_list,mhc_xlsx,output):
    gene_dict = {}
    g_dict ={}
    with open(gene_list,"r")as g:         
        for lines_g in g:
            line_g = lines_g.strip().split("\t")
            chr = line_g[0]
            start = line_g[1]
            end =line_g[2]
            gene =line_g[3]
            index=chr+":"+start+":"+end
            gene_dict[index]=gene
    with open(output+".xls","w") as o:
        with open(output+".txt","w") as ot: 
            with open(mhc_xlsx,"r")as m:
                i = 1         
                for lines in m:
                    if i == 1:
                        i = i + 1
                        o.write(lines)
                        pass
                    elif i ==2:
                        i = i + 1
                        o.write("Pos\tPeptide\tID\tGene\tcore\ticore\tEL-score\tEL_Rank\tAve\tNB")
                        ot.write("Pos\tPeptide\tID\tGene\tcore\ticore\tEL-score\tEL_Rank\tAve\tNB")
                        pass
                    else:
                        line = lines.strip().split("\t")
                        core = line[2]
                        chrid = core.split("_")[0]
                        pos = core.split("_")[1]
                        out_gene = "--"
                        if core not in g_dict:
                            for key in gene_dict.keys():
                                if out_gene != "--":
                                    break
                                else:
                                    chr = key.split(":")[0]
                                    if chr == chrid:
                                        start = key.split(":")[1]
                                        end = key.split(":")[2]
                                        if int(start)<=int(pos)<=int(end):
                                            out_gene = gene_dict[key]
                                            g_dict[core]=out_gene
                        else:
                            out_gene=g_dict[core]
                        out_end = line[3:]
                        out_str = "\t".join(out_end)
                        o.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+out_gene+"\t"+out_str+"\n")
                        ot.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+out_gene+"\t"+out_str+"\n")
    
def main():
    if len(sys.argv) != 4:
        exit("ERROR: This program accepts exactly three arguments: gene_list,mhc_xlsx,output. Exiting...")
    
    print("EXTRACTING FILE " + sys.argv[1])
    
    translate(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
