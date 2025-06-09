#!/usr/bin/env Rscript
library(argparser)
library(tidyr)
library(tidygraph)
library(ggraph)
times <- Sys.time()
argv <- arg_parser("plot_genemania.R")
argv <- add_argument(argv, "--infile", help = "genemania result file")
argv <- add_argument(argv, "--outfile", help = "the output file")
argv <- parse_args(argv)

df <- read.delim(argv$infile, sep = "\t", header = TRUE)
tabs <- table(df$Type)
genes <- sapply(names(tabs),function(x){
    temp <- unlist(df[which(df$Type==x),c("Gene1","Gene2")])
    length(temp[!duplicated(temp)])
})
out <- data.frame(Type=names(tabs),Count=as.numeric(tabs),Genes=genes)
write.table(out,argv$outfile,col.names = T,row.names = F,sep="\t")
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
q()
