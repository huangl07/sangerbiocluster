# -*- coding: utf-8 -*-
# @Last-edit Time 2022/10/12
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com
import re
import argparse

parse = argparse.ArgumentParser(description="生成anno.summary")
parse.add_argument("-i",
                   "--anno",
                   help="输入all_anno_detail.xls文件",
                   required=True)
parse.add_argument("-g", "--gene", help="输入gene_pos.txt文件")
parse.add_argument("-o", "--output", help="输出结果目录", required=True)
args = vars(parse.parse_args())
anno_detail = args["anno"]
gene_pos = args["gene"]
anno_summary = args["output"]
# anno_detail = "all_anno_detail.10000.xls"
# anno_summary = "anno.summary.xls"

gene_info = {}
if gene_pos:
    with open(gene_pos, "r") as f:
        f.readline()
        for line in f:
            item = line.strip().split("\t")
            gene_id = item[0]
            tran_id = item[1]
            chr = item[2]
            start = item[3]
            end = item[4]
            gene_info[gene_id] = chr + ":" + start + ":" + end
            gene_info[tran_id] = chr + ":" + start + ":" + end
with open(anno_detail, "r") as f, open(anno_summary, "w") as w:
    header = "GeneID\tNRID\tNRANNO\tUniID\tUniANNO\tKoID\tKoanno\tGOTERM\tGOANNO\t"
    # header += "EGGNOG\tEGGNOG_ANNO\tPfamID\tPfamAnno\tInterProAccession\tInterProAnno\n"
    header += "EGGNOG\tEGGNOG_ANNO\tPfamID\tPfamAnno\n"
    w.write(header)
    f.readline()
    for line in f:
        item = line.strip().split("\t")
        rna_gene_id = item[0]
        try:
            rna_cog = item[6]
        except:
            rna_cog = ""
        try:
            rna_cog_des = item[7]
        except:
            rna_cog_des = ""
        try:
            rna_ko_paths = item[10]
        except:
            rna_ko_paths = ""
        try:
            rna_pfam = item[11]
        except:
            rna_pfam = ""
        try:
            rna_go = item[12]
        except:
            rna_go = ""
        try:
            rna_nr = item[13]
        except:
            rna_nr = ""
        try:
            rna_swissprot = item[14]
        except:
            rna_swissprot = ""
        anno = []
        tran_id = item[1]
        gene_name = item[3]
        try:
            pos = gene_info[tran_id]
        except:
            try:
                pos = gene_info[rna_gene_id]
            except:
                pos = "--:--:--"
        anno.append(rna_gene_id + ":" + tran_id + ":" + pos)
        nr_id = "--"
        nr_anno = "--"
        p1 = re.compile(r"([^(]*)\((.*)\)$")
        nr_info = re.findall(p1, rna_nr)
        if len(nr_info) > 0:
            nr_id = nr_info[0][0]
            nr_anno = nr_info[0][1]
        anno.append(nr_id)
        anno.append(nr_anno)
        unni_id = "--"
        unni_anno = "--"
        swiss_info = re.findall(p1, rna_swissprot)
        if len(swiss_info) > 0:
            unni_id = swiss_info[0][0]
            unni_anno = swiss_info[0][1]
        anno.append(unni_id)
        anno.append(unni_anno)
        if rna_ko_paths == "":
            anno.append("--")
            anno.append("--")
        else:
            ko_ids, ko_annos = [], []
            ko_info = rna_ko_paths.split("; ")
            for k in ko_info:
                # ko_info_ = re.findall(p1, k)
                # if ko_info_:
                #     map_id = ko_info_[0][0]
                #     ko_id = map_id.replace("map", "ko")
                #     # ko_ids.append(ko_info_[0][0])
                #     ko_ids.append(ko_id)
                #     ko_annos.append(ko_info_[0][1])
                map_id = k.split("(")[0]
                ko_id = map_id.replace("map", "ko")
                ko_ids.append(ko_id)
                p2 = re.compile(r'^[^\(]*\((.*)\)$')
                ko_info_ = re.findall(p2, k)
                if ko_info_:
                    ko_annos.append(ko_info_[0])
            if len(ko_ids) > 0:
                anno.append(",".join(ko_ids))
                anno.append(":".join(ko_annos))
            else:
                anno.append("--")
                anno.append("--")
        if rna_go == "":
            anno.append("--")
            anno.append("--")
        else:
            go_ids, go_annos = [], []
            go_info = rna_go.split("; ")
            for g in go_info:
                g_id = g.split("(")[0]
                go_ids.append(g_id)
                p2 = re.compile(r'^[^\(]*\((.*)\)$')
                go_info_ = re.findall(p2, g)
                if go_info_:
                    g_anno = go_info_[0]
                    g_anno_ = g_anno.split(":",1)
                    if len(g_anno_) > 1:
                        go_annos.append(g_anno_[1])
                    else:
                        go_annos.append(g_anno)
            if len(go_ids) > 0:
                anno.append(",".join(go_ids))
                anno.append(":".join(go_annos))
            else:
                anno.append("--")
                anno.append("--")
        if rna_cog == "":
            anno.append("--")
            anno.append("--")
        else:
            cog, cog_annos = [], []
            cog_info = rna_cog.split("; ")
            cog_info1 = rna_cog_des.split("; ")
            for c in cog_info:
                cog_info_ = re.findall(p1, c)
                if cog_info_:
                    cog.append(cog_info_[0][1])
            for c in cog_info1:
                cog_annos.append(c)
            if len(cog) > 0:
                anno.append(",".join(cog))
                anno.append(":".join(cog_annos))
            else:
                anno.append("--")
                anno.append("--")
        if rna_pfam == "":
            anno.append("--")
            anno.append("--")
        else:
            pfam_ids, pfam_annos = [], []
            pfam_info = rna_pfam.split("; ")
            for p in pfam_info:
                p_id = p.split("(")[0]
                pfam_ids.append(p_id)
                p2 = re.compile(r'^[^\(]*\((.*)\)$')
                pfam_info_ = re.findall(p2, p)
                if pfam_info_:
                    p_anno = pfam_info_[0]
                    p_anno_ = p_anno.split(":",1)
                    if len(p_anno_) > 1:
                        pfam_annos.append(p_anno_[1])
                    else:
                        pfam_annos.append(p_anno)
            if len(pfam_ids) > 0:
                anno.append(",".join(pfam_ids))
                anno.append(":".join(pfam_annos))
            else:
                anno.append("--")
                anno.append("--")
        w.write("\t".join(anno) + "\n")
