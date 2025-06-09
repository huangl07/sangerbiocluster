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
library(AnnotationForge)
library(getopt)

###传参信息
spec <- matrix(c(
  'annofile', 'i', 0, 'character',
  'outpath', 'o', 0, 'character',
  'term2gene', 'n', 0, 'character',
  'term2name', 'm', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--annofile   输入整理好的GO和KEGG注释列表
	--outpath    输出orgdb构建路径
	--term2gene  输出KEGG pathway和对应的transcript_id信息
	--term2name  输出KEGG pathway和对应的通路描述信息
	--help       usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help) )      { print_usage(spec) }
if ( is.null(opt$annofile) )   { print_usage(spec) }
if ( is.null(opt$outpath) )    { print_usage(spec) }
if ( is.null(opt$term2gene) )  { print_usage(spec) }
if ( is.null(opt$term2name) )  { print_usage(spec) }

###导入GO和KEGG的结果
anno_result <- read.table(opt$annofile, header = TRUE, sep = "\t", col.names = c("transcript_id", "KO_id", "ko_id", "ko_anno", "go_term", "go_anno"),
                          colClasses=c("character", "character", "character", "character", "character", "character"), quote = "")

###生成gene_info列表,Gene_name信息就和GID一样即可
gene_name <- anno_result %>% 
  dplyr::select(transcript_id) %>% 
  distinct()
gene_info <- cbind(GID = gene_name, Gene_Name = gene_name)
colnames(gene_info) <- c("GID", "Gene_Name")

###整理GO的结果数据
GO_result <- anno_result %>% 
  dplyr::select(transcript_id, go_term) %>%
  as_tibble() %>% 
  filter(go_term != "", go_term != "--")

###生成gene2go列表
##将GO_result宽数据变长数据，目前是以","符号
gene2go <- GO_result %>% 
  separate_rows(go_term, sep=",",convert = F) %>%
  mutate(EVIDENCE = "IEA")
colnames(gene2go) <- c("GID", "GO", "EVIDENCE")

###生成gene2ko列表

#gene2ko <- anno_result %>% 
#  dplyr::select(transcript_id, KO_id) %>% 
#  filter(KO_id != "--", KO_id != "") %>% 
#  #distinct() %>% 
#  separate_rows(KO_id, sep=",",convert = F) %>% 
#  distinct()
#colnames(gene2ko) <- c("GID", "KO")


###提取GID与Pathway信息(ko值),生成gene2pathway列表
gene2pathway <- anno_result %>% 
  dplyr::select(transcript_id, ko_id) %>% 
  filter(ko_id != "--", ko_id != "") %>% 
  distinct() %>%
  separate_rows(ko_id, sep=",",convert = F)
colnames(gene2pathway) <- c("GID", "pathway")

###建库,在outputDir路径下生成自定义Orgdb包，名字org.Ddemo.eg.db
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               #ko=gene2ko,
               #pathway=gene2pathway,
               maintainer='XY<yuan.xu@majorbio.com>',
               author='yuan.xu',
               version="0.1" ,
               outputDir=opt$outpath, 
               tax_id="1",
               genus="demo",
               species="demo",
               goTable="go")

term2gene <- anno_result %>% 
  dplyr::select(ko_id, transcript_id) %>% 
  filter(ko_id != "--", ko_id != "") %>% 
  distinct() %>%
  separate_rows(ko_id, sep=",",convert = F)
write.table(term2gene, opt$term2gene, sep="\t",col.names = FALSE, row.names = FALSE, quote=FALSE)

term2name <- anno_result %>% 
  dplyr::select(ko_id, ko_anno) %>% 
  filter(ko_id != "--", ko_id != "") %>% 
  mutate(ko_id = gsub(pattern = ",", replacement = ":", .$ko_id,)) %>% 
  separate_rows(ko_id, ko_anno, sep = ":")
write.table(term2name, opt$term2name, sep="\t",col.names = FALSE, row.names = FALSE, quote=FALSE)
