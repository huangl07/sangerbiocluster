#!/usr/bin/env Rscript
times <- Sys.time()
library(getopt)
###传参信息
spec <- matrix(c(
  'seqgz', 's', 0, 'character',
  'outpath', 'o', 0, 'character',
  'samplename', 'n', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--seqgz    输入data.bin50_seqz.gz
	--outpath
    --samplename
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$seqgz))   { print_usage(spec)}
if ( is.null(opt$outpath))  { print_usage(spec) }
if ( is.null(opt$samplename))  { print_usage(spec) }

library(sequenza)
setwd(opt$outpath)
seqz.data <- read.seqz(opt$seqgz)
#Inference of cellularity and ploidy
chromosome_list<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
test <-sequenza.extract(opt$seqgz,chromosome.list=chromosome_list)
CP.example <- sequenza.fit(test)
#Results of model fitting
sequenza.results(sequenza.extract = test, cp.table = CP.example,sample.id = opt$samplename, out.dir=opt$outpath)
times <- Sys.time() 