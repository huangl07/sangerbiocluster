#!/usr/bin/env Rscript

suppressMessages({
library(tidyverse)
library(argparser)})
argv <- arg_parser('处理final_anno文件')
argv <- add_argument(argv,"--infile", help="the input file")
argv <- add_argument(argv,"--outfile", help="the output file")
argv <- parse_args(argv)
infile <- argv$infile
outfile <- argv$outfile

data <- read.table(infile,sep="\t",head=T,quote="",comment.char="")
out_data <- data %>%
  #select(-InterProAccession, -InterProAnno, -KOID) %>%
  mutate(Rank=case_when(Rank != "--" ~ paste0("(",Rank,")"),
                        Rank == "--" ~ "--"))

write.table(file=outfile, out_data, sep="\t",quote=F,row.names=F)
