library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'out','o',1,'character',
	'pop','p',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript emap.R --mark --out --pop
	
Usage:
	--mark	map file
	--out	out dir
	--pop	pop type
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }
if(grepl("RIL",opt$pop) | grepl("Ril",opt$pop) | grepl("Ri",opt$pop)){opt$pop="riself"}
if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}
library('qtl');



opt$pop=tolower(opt$pop)
if(opt$pop == "cp" | opt$pop == "CP"){
	map = paste(opt$mark,"sexAver.map",sep=".")
	if(!file.exists(map)){
		map=paste(opt$mark,"map",sep=".")
	}
	d<-read.cross(genfile=paste(opt$mark,"loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=map,format="mapqtl",crosstype="4way")
	setwd(opt$out)
	d<-jittermap(d)
	d<-est.rf(d)
	chrname<-chrnames(d);
	for(i in c(1:length(chrname))){
		pdf(paste("chr",chrname[i],".heatMap.sexAver.pdf",sep=""))
		plotRF(d,chr=i,lmax=50,"both")
		dev.off()
		png(paste("chr",chrname[i],".heatMap.sexAver.png",sep=""))
		plotRF(d,chr=i,lmax=50,"both")
		dev.off()
	}
	
}else{
	library('ASMap');
	d<-read.cross(file=opt$mark,format="csvr",crosstype=opt$pop)
	setwd(opt$out)
	d<-jittermap(d)
	d<-est.rf(d)

	chrname<-chrnames(d);
	for(i in c(1:length(chrname))){
		pdf(paste("chr",i,".heatMap.pdf",sep=""))
		plotRF(d,chr=i,what='both')
		dev.off()
		png(paste("chr",i,".heatMap.png",sep=""))
		plotRF(d,chr=i,what='both')
		dev.off()
	}
}

escaptime=Sys.time()-times;
print(escaptime)
