#!/usr/bin/env Rscript
start_time <- Sys.time()
library(getopt)
library(ConsensusClusterPlus)
options(bitmapType = "cairo")
spec <- matrix(c(
  'maf', 'm', 0, 'character',
  'out', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--maf    输入maf表格
    --out    输出基因突变全景图
	--help   usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$maf))   { print_usage(spec) }
if ( is.null(opt$out))   { print_usage(spec) }
library(maftools)
laml<-read.maf(opt$maf)
infer<-inferHeterogeneity(maf=laml,vafCol="i_TumorVAF_WU")
pdf(paste0(opt$out,".pdf"),height=10,width=10)
pv <- plotClusters(clusters=infer)
dev.off()
png(paste0(opt$out,".png"),height=1500,width=1500,unit="px",res=300)
pv <- plotClusters(clusters=infer)
dev.off()
x = tmb(maf=laml)
write.table (x, file =paste0(opt$out,"TMB.txt"), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
write.table (infer$clusterMeans, file =paste0(opt$out,"cluster_means.txt"), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
write.table (infer$clusterData, file =paste0(opt$out,"cluster_data.txt"), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
