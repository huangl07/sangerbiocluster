#! /mnt/ilustre/app/pub/R/bin/Rscript
library(getopt)
options(bitmapType='cairo')
opt = getopt(matrix(c(
'groupid','g',1,'character',
'out','o',2,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
cat("This script is used to plot genetic map
Usage	Rscript map-plot.R for fig 3-9[options]
Options:
	--groupid, -g	trt file, forced 
	--out, -o 	out file, forced
	-h, --help	print display this help and exit
")
q(status=1);
}
times<-Sys.time()
if(is.null(opt$groupid)){usage()}
if(is.null(opt$out)){usage()}


library(maftools)

laml.gistic = readGistic(gisticAllLesionsFile = paste(opt$groupid,"all_lesions.conf_75.txt",sep="."), gisticAmpGenesFile = paste(opt$groupid,"amp_genes.conf_75.txt",sep="."), gisticDelGenesFile = paste(opt$groupid,"del_genes.conf_75.txt",sep="."), gisticScoresFile = paste(opt$groupid,"scores.gistic",sep="."), isTCGA = FALSE)

if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out)
 pdf(paste(opt$groupid,"gisticChrom.pdf",sep="."),height=9,width=16);
 gisticChromPlot(gistic = laml.gistic, markBands = "all",cytobandOffset=0.25);
 dev.off()
 png(paste(opt$groupid,"gisticChrom.png",sep="."),height=900,width=1600);
 gisticChromPlot(gistic = laml.gistic, markBands = "all",cytobandOffset=0.25);
 dev.off()


 pdf(paste(opt$groupid,"gistic.oncoplot.pdf",sep="."),height=9,width=16);
 gisticOncoPlot(gistic = laml.gistic, showTumorSampleBarcodes=T,top=20);
 dev.off()
 png(paste(opt$groupid,"gistic.oncoplot.png",sep="."),height=900,width=1600);
 gisticOncoPlot(gistic = laml.gistic,showTumorSampleBarcodes=T,top=20);
 dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
