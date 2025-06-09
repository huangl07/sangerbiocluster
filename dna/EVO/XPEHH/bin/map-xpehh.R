library('getopt');
times<-Sys.time();

spec = matrix(c(
	'map','v',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--map	the input hapmap file
	--output	the output dir
	--help		usage
\n")
	q(status=1);
}
map<-read.table(opt$map,sep="\t",head =F)
map$V3=map$V4/max(map$V4) * 100
write.table(map,file=opt$output,col.names=F,row.names=F,quote=F)



escaptime=Sys.time()-times
print("Done!")
print(escaptime)
