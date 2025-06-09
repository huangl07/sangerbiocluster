#!/usr/bin/env Rscript
###安装AnnotationForge包
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("AnnotationForge", quietly = TRUE)){
  BiocManager::install("AnnotationForge")
}
if (!requireNamespace("GO.db", quietly = TRUE)){
  BiocManager::install("GO.db")
}

###加载R包
library(tidyverse)
library(purrr)
library(AnnotationForge)
library(getopt)

###传参信息
spec <- matrix(c(
  'gofile', 'r', 0, 'character',
  'keggfile', 'g', 0, 'character',
  'outpath', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--gofile    输入GO结果列表
	--keggfile   输入处理完规整的kegg结果列表
	--outpath   输出orgdb构建路径
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$gofile))   { print_usage(spec)}
if ( is.null(opt$keggfile))  { print_usage(spec) }
if ( is.null(opt$outpath))  { print_usage(spec) }

###导入GO的结果数据
GO_result <- read.delim(opt$gofile, header = FALSE, col.names = c("GID", "GO")) %>% 
  as_tibble() %>%
  na.omit()

###导入KEGG的结果数据
KEGG_result <- read.delim(opt$keggfile, header = TRUE, col.names = c("GID", "KO", "pathway_description", "pathway")) %>% as_tibble()

###生成gene_info列表,Gene_name信息就和GID一样即可
gene_name <- GO_result %>% 
  dplyr::select(GID) %>% 
  distinct()
gene_info <- cbind(GID = gene_name, Gene_Name = gene_name)
colnames(gene_info) <- c("GID", "Gene_Name")

###生成gene2go列表
##将GO_result宽数据变长数据，目前是以"; "符号
# 把一行GO结果变成多行(有更好的写法)
longer_go_result<- function(x,y){
  y_split <- str_split(y, "; ")[[1]]
  len_y <- length(y_split)
  x_rep <- rep(x, len_y)
  GO_tibble <- tibble(GID = x_rep, GO = y_split)
  return(GO_tibble)
}
#gene2go <- map2_dfr(GO_result$GID, GO_result$GO, longer_go_result)
gene2go <- GO_result %>% 
  separate_rows(GO, sep="; ",convert = F)
##去除"--"和""的结果，保留GO结果;加上
gene2go <- gene2go %>% 
  filter(startsWith(GO,"GO")) %>%
  mutate(EVIDENCE = "IEA")

###生成gene2ko列表
gene2ko <- KEGG_result %>% 
  dplyr::select(GID, KO) %>% 
  filter(KO != "-") %>% 
  distinct()

###提取GID与Pathway信息(ko值),生成gene2pathway列表
gene2pathway <- KEGG_result %>% 
  dplyr::select(GID, pathway) %>% 
  filter(pathway != "-") 

###建库,在outputDir路径下生成自定义Orgdb包，名字org.Ddemo.eg.db
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               maintainer='XY<yuan.xu@majorbio.com>',
               author='yuan.xu',
               version="0.1" ,
               outputDir=opt$outpath, 
               tax_id="1",
               genus="demo",
               species="demo",
               goTable="go")

