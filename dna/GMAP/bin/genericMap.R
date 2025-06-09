times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'input','a',0,'character',
	'output','b',0,'character',
    'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript QC.r --base  --qual  --key  --od
	
Usage:
	--input	    base distribution file
	--output	base quality file
    --help		usage
\n")
	q(status=1);
}

if (!is.null(opt$help) ) { print_usage(spec) }
if(is.null(opt$input)){print_usage(spec)}
if(is.null(opt$output)){print_usage(spec)}

a<-read.table(opt$input,sep="\t",skip=1)
myfun<-function(l)(if(max(l,na.rm=T) > 200 || max(l,na.rm=T) < 80){l/max(l,na.rm=T)* runif(1,100,150)}else{l})
#sapply(a[-1],myfun)
result=data.frame(a[1],sapply(a[-1],myfun))
write.table(file=opt$output,result,sep="\t",quote=F,row.names=F,col.names=F)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
