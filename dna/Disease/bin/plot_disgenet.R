#!/usr/bin/env Rscript
library(argparser)
library(tidyr)
library(tidygraph)
library(ggraph)
times <- Sys.time()
argv <- arg_parser("plot_disgenet.R")
argv <- add_argument(argv, "--genefile", help = "disgenet gene result file")
argv <- add_argument(argv, "--varfile", help = "disgenet variant result file")
argv <- add_argument(argv, "--outfile", help = "the output file")
argv <- parse_args(argv)

gene_df <- read.delim(argv$genefile, sep = "\t", header = TRUE)
var_df <- read.delim(argv$varfile, sep = "\t", header = TRUE)
if(nrow(gene_df)+nrow(var_df)>0){
    disease_df <- rbind(gene_df[, c("diseaseId", "diseaseName")],var_df[, c("diseaseId", "diseaseName")])
    disease_df <- disease_df[!duplicated(disease_df), ]
    nodes <- data.frame(
        name = c(as.character(gene_df$geneId), disease_df$diseaseId, var_df$variantId),
        symbol = c(gene_df$geneName, disease_df$diseaseName, var_df$variantId),
        nodeType = c(rep_len("gene", nrow(gene_df)), rep_len("disease", nrow(disease_df)), rep_len("SNV", nrow(var_df)))
    )
    nodes <- nodes[!duplicated(nodes),]
    edges <- data.frame(
        from = c(as.character(gene_df$geneId), var_df$variantId),
        to = c(gene_df$diseaseId, var_df$diseaseId),
        link = c(gene_df$score, var_df$score)
    )
    edges <- edges[!duplicated(edges),]
    gr <- tbl_graph(edges = edges, nodes = nodes)
    plot <- ggraph(gr, "graphopt") +
        geom_edge_link(aes(edge_width = link), color = "grey50") +
        geom_node_point(aes(color = nodeType), size = 25) +
        geom_node_text(aes(label = symbol)) + coord_cartesian(clip = "off") +
        scale_color_manual(name = "Node", breaks = c("gene", "disease","SNV"), values = c("orange2", "chartreuse3","mediumpurple2")) +
        coord_cartesian(clip = "off") +
        theme(
            panel.border = element_blank(), panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_blank(), axis.ticks = element_blank()
        )
    ggsave(paste0(argv$outfile, ".pdf"), plot, dpi = 600, width = 10, height = 10)
    ggsave(paste0(argv$outfile, ".png"), plot, dpi = 600, width = 10, height = 10)
} else {
    print("No disgenet record found!")
}
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
q()
