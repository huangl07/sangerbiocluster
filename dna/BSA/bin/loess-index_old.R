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
data<-read.table(opt$infile,head=TRUE,comment.char="^")
chrlist=levels(as.factor(data$X.chr))
loess0<-function(chr){
   df1 =data[data$X.chr == chr,]
   if(nrow(df1) < 100){
        return=data.frame(chr=df1$X.chr,pos=df1$pos,loess=1)
   }else{
        result=loess(delta ~ pos,data=df1) %>% predict(df1$pos)
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
write.table(file=opt$outfile,result,sep="\t",row.names=F)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
