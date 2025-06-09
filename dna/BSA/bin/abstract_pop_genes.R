#!/usr/bin/env Rscript
times <- Sys.time()
if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)
###传参信息
spec <- matrix(c(
  'regionfile', 'r', 0, 'character',
  'genefile', 'g', 0, 'character',
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
	--outfile   输出文件
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$regionfile))   { print_usage(spec)}
if ( is.null(opt$genefile))  { print_usage(spec) }
if ( is.null(opt$outfile))  { print_usage(spec) }

pacman::p_load(tidyverse, purrr)

pop_region<-read.table(opt$regionfile)
colnames(pop_region)=c("chr","pos_start","pos_end","nGene")
###导入所有的gene_list
all_gene_list <- read.delim(opt$genefile, header = TRUE, col.names = c("chr", "transcript_start", "transcript_end", "transcript_id",
                                                                       "gene_start", "gene_end", "gene_id"),
                            colClasses=c("character", "numeric", "numeric", "character", "numeric", "numeric", "character")) %>%
  as_tibble()

abstract_df <- function(x, pos_start, pos_end, all_gene_list){
  df_abstract <- all_gene_list %>%
    filter(chr == x, gene_end > pos_start, gene_start < pos_end) %>%
    mutate(REGION = paste(x, pos_start, pos_end, sep=":")) %>%
    select(REGION, chr, gene_start, gene_end, gene_id, transcript_id)
  return(df_abstract)
}
df_abstract_all <- pmap_dfr(pop_region[1:3], ~abstract_df(..1, ..2, ..3, all_gene_list))


###结果输出
write.table(df_abstract_all, opt$outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
