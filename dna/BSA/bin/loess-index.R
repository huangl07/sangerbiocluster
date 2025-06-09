#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
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
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) )   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outfile))  { print_usage(spec) }
opt$winsize=as.numeric(opt$winsize)
opt$power=as.numeric(opt$power)
opt$stepsize=as.numeric(opt$step)
library(dplyr)
library(caret)
	ctrl=trainControl(method="cv",number=5)
	grid <- expand.grid(span = seq(0.5, 0.9, len = 10), degree = 1)
data<-read.table(opt$infile,head=TRUE,comment.char="^")
chrlist=levels(as.factor(data$X.chr))


loess0<-function(chr){
	print(chr)
	df1 =data[data$X.chr == chr,]
	if(nrow(df1) < 50){
        return=data.frame(chr=df1$X.chr,pos=df1$pos,loess=0)
	}else{
		model <- train(delta ~ pos, data = df1, method = "gamLoess", tuneGrid=grid, trControl = ctrl)
        result=loess(delta ~ pos,data=df1,span=model$bestTune$span) %>% predict(df1$pos)
        return=data.frame(chr=df1$X.chr,pos=df1$pos,loess=result)
	}
	return
}
rlist=lapply(chrlist,loess0)


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
