#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'winsize','w',0,'character',
	'stepsize','s',0,'character',
	'abs','a',0,'character',
	'method','m',0,'character',
	'power','p',0,'character',
	'minmarker','x',0,'numeric',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--infile	the input hapmap file
	--outfile	the trait file
	--winsize	the window size
	--stepsize	the window sliding step
	--method	the sliding method default bp or num
	--power	the power of ed calc
	--abs	the delta is abs or not
	--minmarker	min marker number in a sliding window
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) )   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outfile))  { print_usage(spec) }
if ( is.null(opt$winsize))  {opt$winsize=1000000}
if ( is.null(opt$stepsize))  {opt$stepsize=10000}
if (is.null(opt$method)){opt$method="bp"}
if(is.null(opt$power)){opt$power=4}
if(is.null(opt$minmarker)){opt$minmarker=10}
opt$winsize=as.numeric(opt$winsize)
opt$power=as.numeric(opt$power)
opt$stepsize=as.numeric(opt$stepsize)
print(opt$stepsize)
library(dplyr)
library("Rcpp")
argv <- commandArgs(trailingOnly=F)
temp <- '--file='
script.abspath <- sub(temp, "", argv[grep(temp, argv)])
script.dirname <- dirname(script.abspath)
print(script.abspath)
sourceCpp(file.path(script.dirname,"countSNPs-step.cpp"))
sourceCpp(file.path(script.dirname,"countDepth-step.cpp"))

data<-read.table(opt$infile,head=TRUE,comment.char="^")
if(ncol(data) == 4){
	if(opt$method == "num"){
		data<-data %>%
			group_by(X.chr) %>%
			mutate(psites = seq_along(pos) )
	}else{
		data$psites=data$pos
	}
	data$pos1=ceiling(data$psites/opt$stepsize)*opt$stepsize-floor(opt$winsize/2)
	data$pos2=ceiling(data$psites/opt$stepsize)*opt$stepsize+floor(opt$winsize/2)
	data$pos1[data$pos1 <=0]=ifelse(opt$method == "num", 1, 0)
	# data$pos2[data$pos1 <=0]=opt$winsize
	df1=data[data$psites-data$pos1 < opt$stepsize & data$pos1 != ifelse(opt$method == "num", 1, 0),]
	df1$pos1=df1$pos1-opt$stepsize
	df1$pos2=df1$pos2-opt$stepsize
	data=rbind(data,df1)

	data<-data %>% group_by(X.chr) %>%
		mutate(mdepth=mdepths(POS = psites, DEPTH=depth,pos1=pos1,pos2=pos2))
	data<-data %>% group_by(X.chr) %>%
		mutate(ncount=countSNPs_cpp(POS = psites,pos1=pos1,pos2=pos2))
	data<-data %>% group_by(X.chr) %>%
		mutate(slidingD=mdepths(POS = psites, DEPTH=delta,pos1=pos1,pos2=pos2))
	data$mdepth[data$ncount < opt$minmarker]=mean(data$mdepth)
	data$slidingD[data$ncount < opt$minmarker] =mean(data$slidingD)
	data$mdepth=ceiling(data$mdepth)
	df=data.frame(X.chr=data$X.chr,pos1=data$pos1,pos2=data$pos2,mdepth=data$mdepth,slidingD=data$slidingD)
	df=df[!duplicated(df),]
	write.table(file=paste(opt$outfile,"sliding.result",sep="."),df,quote=F,row.names=F,sep="\t")
	write.table(file=paste(opt$outfile,"sliding.detail",sep="."),data,quote=F,row.names=F,sep="\t")
}else{
	data$depth=data$n1+data$n2+data$n3+data$n4
	data$ref_depth = data$n1 + data$n3
	data$alt_depth = data$n2 + data$n4
	data$ED=sqrt((data$n1/data$ref_depth-data$n2/data$alt_depth)^2+(data$n3/data$ref_depth-data$n4/data$alt_depth)^2)
	data$ED4=data$ED^opt$power
	if(opt$method == "num"){
		data<-data %>%
			group_by(X.chr) %>%
			mutate(psites = seq_along(pos) )
	}else{
		data$psites=data$pos
	}
	data$pos1=ceiling(data$psites/opt$stepsize)*opt$stepsize-floor(opt$winsize/2)
	data$pos2=ceiling(data$psites/opt$stepsize)*opt$stepsize+floor(opt$winsize/2)
	data$pos1[data$pos1 <=0]=data$pos1[data$pos1 <=0]=ifelse(opt$method == "num", 1, 0)
	# data$pos2[data$pos1 <=0]=opt$winsize
	df1=data[data$psites-data$pos1 < opt$stepsize & data$pos1 !=ifelse(opt$method == "num", 1, 0),]
	df1$pos1=df1$pos1-opt$stepsize
	df1$pos2=df1$pos2-opt$stepsize
	data=rbind(data,df1)
	print(max(data$pos1))
	data<-data %>% group_by(X.chr) %>%
		mutate(slidingI1=mdepths(POS = psites, DEPTH=index1,pos1=pos1,pos2=pos2))
	data<-data %>% group_by(X.chr) %>%
		mutate(slidingI2=mdepths(POS = psites, DEPTH=index2,pos1=pos1,pos2=pos2))
	data<-data %>% group_by(X.chr) %>%
		mutate(mdepth=mdepths(POS = psites, DEPTH=depth,pos1=pos1,pos2=pos2))
	data <- data %>% group_by(X.chr) %>%
		mutate(slidingED=mdepths(POS=psites,DEPTH=ED4,pos1=pos1,pos2=pos2))
	data<-data %>% group_by(X.chr) %>%
		mutate(ncount=countSNPs_cpp(POS = psites,pos1=pos1,pos2=pos2))
  data<-data %>% group_by(X.chr) %>%
  mutate(slidingD=mdepths(POS = psites, DEPTH=delta,pos1=pos1,pos2=pos2))
	if (is.null(opt$abs)){
	data$slidingD2=data$slidingI1-data$slidingI2
	}else{
	data$slidingD=abs(data$slidingD)
	}
	data$mdepth[data$ncount < opt$minmarker]=mean(data$mdepth)
	data$slidingI2[data$ncount < opt$minmarker] =mean(data$slidingI2)
	data$slidingI1[data$ncount < opt$minmarker] =mean(data$slidingI1)
	data$slidingD[data$ncount < opt$minmarker] =mean(data$slidingD)
	data$mdepth=ceiling(data$mdepth/2)
	df=data.frame(X.chr=data$X.chr,pos1=data$pos1,pos2=data$pos2,slidingI1=data$slidingI1,slidingI2=data$slidingI2,mdepth=data$mdepth,slidingD=data$slidingD,slidingED=data$slidingED)
	df$pos=(df$pos1+data$pos2)/2
	df=df[!duplicated(df),]
	write.table(file=paste(opt$outfile,"sliding.result",sep="."),df,quote=F,row.names=F,sep="\t")
	write.table(file=paste(opt$outfile,"sliding.detail",sep="."),data,quote=F,row.names=F,sep="\t")
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
