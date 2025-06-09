#!/usr/bin/env Rscript
times <- Sys.time()

library(getopt)
###传参信息
spec <- matrix(c(
  'infile', 'i', 0, 'character',
  'ccol', 'c', 0, 'character',
  'pos','a',0,'character',
  'loess','l',0,'character',
  'CI','x',0,'character',
  "quantile","q",0,'character',
  'outfile', 'o', 0, 'character',
  'number', 'p', 0, 'integer',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--infile    the input pop.bootstrap.result file
	--ccol   the colnames for chromsome
  --pos    the colnames for position
  --loess  the colnames for loess
  --CI     the colnames for threshold
  --quantile  the quantile threshold for ED
	--outfile   the abstract_result file
	--number   number of continuous default 10
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outfile))  { print_usage(spec) }
if ( is.null(opt$threshold))  { opt$threshold=10 }

if(is.null(opt$ccol)){opt$ccol = "chr"}

if(is.null(opt$pos)){opt$pos = "chr"}
if(is.null(opt$loess)){opt$loess = "chr"}
if(is.null(opt$CI)){opt$CI = "chr"}

if (!require("pacman")){
  install.packages("pacman")
}
library(pacman)
p_load(tidyverse, data.table, purrr)

###读取pop.bootstrap.result文件和指定列名文件
pop_bootstrap_result <- fread(opt$infile) %>% 
  tibble() 
  if(!is.null(opt$quantile)){
      opt$CI="CI"
      pop_bootstrap_result$CI=quantile(pop_bootstrap_result[[opt$loess]],as.numeric(opt$quantile),na.rm=T)
  }
setnames(pop_bootstrap_result, c(opt$ccol,opt$pos,opt$loess,opt$CI), c("chr", "pos", "loess", "CI")) #重命名

###比较loess和CI列,分组并提取子集
tmp <- pop_bootstrap_result$loess > pop_bootstrap_result$CI
pop_bootstrap_result$tmp <- tmp
list_by_chr <- pop_bootstrap_result %>% group_split(chr)
##间隔分组
group_intervals <- function(x, threshold){
  x$ID <- rleid(x$tmp)
  x <- x %>% filter(tmp == TRUE) %>% group_by(ID) %>% filter(n() > threshold) 
  return(x)
}
rle_result <- lapply(list_by_chr,group_intervals, opt$threshold) %>% do.call(rbind, .) 
##输出符合条件的子集
split_result <- rle_result %>% group_by(chr, ID) %>% group_split()
##统计每一组的结果
result_split <- function(x){
  group_start <- 1
  group_end <- dim(x)[1]
  chr <- x$chr[1]
  pos_start <- x$pos[1]
  pos_end <- x$pos[group_end]
  num_snv <- group_end
  return(c(chr = chr, pos_start = pos_start, pos_end = pos_end, num_snv = num_snv))
}
result <- lapply(split_result,result_split) %>% do.call(rbind, .) %>% as_tibble()
##输出最终结果
write.table(result, file = opt$outfile, sep = "\t", quote=FALSE, row.names=FALSE)




escaptime <- Sys.time()-times
print("Done!")
print(escaptime)
q()
