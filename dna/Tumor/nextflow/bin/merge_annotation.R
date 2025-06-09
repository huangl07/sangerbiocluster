#!/usr/bin/env Rscript

if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)
###传参信息
spec <- matrix(c(
  'regionfile', 'r', 0, 'character',
  'genefile', 'l', 0, 'character',
  'gofile', 'g', 0, 'character',
  'keggfile', 'k', 0, 'character',
  'outfile', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--regionfile    输入pop.ED.region列表
	--genefile   输入all_gene_list列表
  --gofile    输入GO结果文件
  --keggfile  输入整理过的KEGG文件
	--outfile   输出文件
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$regionfile))   { print_usage(spec)}
if ( is.null(opt$genefile))  { print_usage(spec) }
if ( is.null(opt$gofile))   { print_usage(spec)}
if ( is.null(opt$keggfile))  { print_usage(spec) }
if ( is.null(opt$outfile))  { print_usage(spec) }

pacman::p_load(tidyverse, purrr, data.table)

abstract_df <- function(x, pos_start, pos_end, all_gene_list){
  df_abstract <- all_gene_list %>%
    filter(chr == x, gene_start < pos_end, gene_end > pos_start) %>%
    mutate(region = str_c(x," ",pos_start," ",pos_end))
  return(df_abstract)
}

fill_na <- function(x){x[is.na(x )]<- "-";x}

###导入pop.ED.region的结果
pop_region <- read.delim(opt$regionfile, header = TRUE, colClasses=c("chr"= "character")) %>% as_tibble()
###导入所有的gene_list
all_gene_list <- read.delim(opt$genefile, header = TRUE, col.names = c("REGION","chr", "gene_start", "gene_end", "gene_id"), colClasses=c("character","character", "numeric", "numeric", "character")) %>% as_tibble()
###导入GO的结果；去除"--"和""的结果，保留GO结果
GO_result <- read.delim(opt$gofile, header = FALSE, col.names = c("GID", "GO")) %>% 
  as_tibble() %>% 
  na.omit() %>% 
  filter(startsWith(GO,"GO"))
###导入KEGG结果，并合并pathway列
KEGG_result <- read.delim(opt$keggfile, header = TRUE, col.names = c("GID", "KO", "pathway_description", "pathway")) %>% 
  as_tibble() %>% 
  group_by(GID, KO) %>% 
  summarise(pathway = str_c(pathway, collapse = "; "), .groups ="drop")
###提取region中的所有基因
df_abstract_all <- pmap_dfr(pop_region[1:3], ~abstract_df(..1, ..2, ..3, all_gene_list))

###连接df_abstract_all和GO结果
df_result <- merge(df_abstract_all, GO_result, all.x=TRUE, by.x = "gene_id", by.y = "GID")

###连接df_abstract_all和KEGG结果
df_result <- merge(df_result, KEGG_result,all.x=TRUE, by.x = "gene_id", by.y = "GID")

###填充缺失值为"-"
df_result <- fill_na(df_result)

###整理列
df_result <- df_result %>% 
  select(region,everything()) %>% 
  relocate(gene_id, .after = gene_end)

###结果输出
write.table(df_result, opt$outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
