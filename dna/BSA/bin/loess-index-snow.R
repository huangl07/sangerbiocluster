#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'mutmap','m',0,'logical',
	'help','h',0,'logical',
    'minmarker','n',1,'numeric'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--infile	the input hapmap file
	--outfile	the trait file
	--mutmap    whether is mutmap
	--help		usage
    --minmarker the minimum marker number (default 5, must > 1)
\n")
	q(status=1);
}
if ( !is.null(opt$help) )   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outfile))  { print_usage(spec) }
if ( is.null(opt$mutmap))   { opt$mutmap <- FALSE }
if ( is.null(opt$minmarker))   { opt$minmarker <- 5 }
opt$winsize=as.numeric(opt$winsize)
opt$power=as.numeric(opt$power)
opt$stepsize=as.numeric(opt$step)
library(dplyr)
library(caret)
ctrl=trainControl(method="cv",number=5)
grid <- expand.grid(span = seq(0.5, 0.9, len = 10), degree = 1)
data<-read.table(opt$infile,head=TRUE,comment.char="^")
chrm <- data%>%group_by(X.chr)%>%summarize(n=n())
data <- data[data$X.chr %in% chrm$X.chr[chrm$n>=opt$minmarker],]
chrlist <- levels(as.factor(data$X.chr))

library(snowfall)
sfInit(parallel=TRUE,cpus=8,slaveOutfile="log.txt")
sfLibrary(caret)
sfLibrary(dplyr)
sfExport("data")
sfExport("chrlist")
sfExport("ctrl")
sfExport("grid")
sfExport("opt")
loess0<-function(chr){
	message(chr)
	df1 =data[data$X.chr == chr,]
	if(nrow(df1) < 100 && !opt$mutmap){
        rdf=data.frame(chr=df1$X.chr,pos=df1$pos,loess=rep(0,nrow(df1)))
	}else{
		if(nrow(df1) > 50000){
			df2 = df1[createDataPartition(df1$pos,p=50000/nrow(df1))$Resample1,]
		}else{
			df2 = df1
		}
		model <- train(delta ~ pos, data = df2, method = "gamLoess", tuneGrid=grid, trControl = ctrl,allow.parallel=T)
        result=predict(model$finalModel,newdata=df1)
        rdf=data.frame(chr=df1$X.chr,pos=df1$pos,loess=result)
	}
	rdf
}
print("haha")
rlist=sfLapply(chrlist,loess0)
sfStop()

names(rlist)=chrlist
df=do.call(rbind,rlist)
colnames(data)[1]="chr"
print("haha")
if("n1" %in% colnames(data)){
	data$mdepth=(data$n1+data$n2+data$n3+data$n4)/2
}else{
	data$mdepth=data$depth
}
result=left_join(data,df,by=c("chr","pos"))
print(colnames(result))
colnames(result)[1]="X.chr"
result$pos1=result$pos
result$pos2=result$pos
result$slidingI1=result$index1
result$slidingI2=result$index2
result$slidingD=result$loess
write.table(file=opt$outfile,result,sep="\t",row.names=F,quote=F)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
