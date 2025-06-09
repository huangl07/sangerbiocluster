#! /mnt/ilustre/app/pub/R/bin/Rscript
library(getopt)
options(bitmapType='cairo')
opt = getopt(matrix(c(
'trt','i',1,'character',
'out','o',2,'character',
'order','d',2,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
cat("This script is used to plot genetic map
Usage	Rscript map-plot.R for fig 3-9[options]
Options:
	-m, --trt	trt file, forced 
	-o, --out 	out file, forced
	-d, --order input order file
	-h, --help	print display this help and exit
")
q(status=1);
}
times<-Sys.time()

if (!is.null(opt$help) ) { usage() }
if (is.null(opt$trt) ) { usage() }
library(dplyr)
library(reshape2)
txt<-read.table(opt$trt,sep="\t",head=T)
if(!is.null(opt$order)){
	order=read.table(opt$order,head=F)
	rownames(txt)=as.character(txt[,1])
	txt=txt[as.character(order$V1),]
	print(head(txt))
}
colnames(txt)[1]="sampleID"
id=colnames(txt)[1];
txt[is.na(txt)]="*"
df=melt(txt,id.vars=c("sampleID"))

print(head(df))
stat=df %>% group_by(variable) %>% summarise(nlevels=nlevels(as.factor(value)))
qname=as.factor(stat$variable[stat$nlevels > 2])
bname=as.factor(stat$variable[stat$nlevels == 2])
bname=as.character(bname)
nd=txt[,c(id,bname)]

myfun<-function(name){as.numeric(as.factor(nd[[name]]))-1}
nnd=as.data.frame(lapply(bname,myfun),col.names=bname)
if(nrow(nnd)!=0){nnd[[id]]=txt[[id]]}

print(qname)
print(bname)
write.table(txt[,c(id,as.character(qname))],file=paste(opt$out,"qtl.txt",sep="."),quote=FALSE,sep=" ",row.name=F)
if(nrow(nnd)!=0){
write.table(nnd[,c(id,bname)],file=paste(opt$out,"btl.txt",sep="."),quote=FALSE,sep=" ",row.name=F)
}else{
	cat("",file=paste(opt$out,"btl.txt",sep="."))
}


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
