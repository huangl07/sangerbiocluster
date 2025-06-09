library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
    'map','m',1,'character',
    'loc','l',1,'character',
    'trt','t',1,'character',
    'out','o',1,'character',
    'num','n',1,'character',
    'pvalue','q',1,'character',
    'popt','p',1,'character',
    'method','e',1,'character',
    'help','h',0,'logical'
    ), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
Usage example: 
    Rscript Rqtl_CP.r --map --loc --trt --out --num
    
Usage:
    --map    map file
    --loc    loc file
    --trt    trt file
    --out    out dir
    --num    pm number
    --popt  CP
    --pvalue    select pvalue
    --method    
    --help        usage
\n")
    q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$popt) ) { opt$popt = "CP"}
if(opt$popt == "CP"){
    if(is.null(opt$map)){print_usage(spec)}
    opt$method="scanone"
}else{
    if (is.null(opt$method)){opt$method="cim"}
}

if ( is.null(opt$map) & opt$popt == "CP") { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }
if ( is.null(opt$out) ) { opt$out="./";}
if ( is.null(opt$pvalue) ) { opt$pvalue=0.05}else{opt$pvalue=as.numeric(opt$pvalue)}
opt$num=as.numeric(opt$num)
if(opt$popt == "CP"){
    d<-read.cross(mapfile=opt$map,genfile=opt$loc,phefile=opt$trt,format="mapqtl",crosstype="4way")
}else{
    d<-read.cross(format="csvsr", genfile=opt$loc, phefile=opt$trt,crosstype=opt$popt)
}
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno)[2]
print(phe.name)
pdf(paste(opt$out,"pheno.pdf",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
png(paste(opt$out,"pheno.png",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
model="normal"
if(opt$method == "scanone"){
scan<-scanone(d,pheno.col=phe.name,model=model);
scan.pm<-scanone(d,pheno.col=phe.name,n.perm=opt$num,model=model);
}else{
    opt$method == "cim"
}
markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)

pm.result<-summary(scan.pm,alpha=opt$pvalue)
detail=data.frame(marker=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod,pm1=pm.result[1,1])
scan.result<-summary(scan,format="tabByCol",threshold=pm.result,drop=1)

if(nrow(scan.result$lod) ==0){scan.result<-summary(scan,format="tabByCol",threshold=max(scan$lod)-0.001,drop=1);pm.result=c(3)}
print(scan.result$lod)
markerid1<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.high)
markerid2<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.low)
qtl<-makeqtl(d,chr=scan.result$lod$chr,pos=scan.result$lod$pos)
if(length(markerid1) > 1){
qtl<-refineqtl(cross=d,qtl=qtl,pheno.col=phe.name)
}
fitqtl<-fitqtl(cross=d,qtl=qtl,pheno.col=phe.name,get.est=TRUE)
var=fitqtl$result.full["Model","%var"]
if(length(qtl$name)  > 1){var=fitqtl$result.drop[,"%var"]}
result=data.frame(scan.result$lod,pos1=markerid1,pos2=markerid2,var=var)
write.table(detail,file=paste(opt$out,"detail.result",sep="."),row.names=F,quote=F)
write.table(result,file=paste(opt$out,"qtl-result.result",sep="."),row.names=F,quote=F)
write.table(file=paste(phe.name,".pm.csv",sep=""),sep="\t",scan.pm,quote=F,row.names=T);

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
