library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'exposure','e',1,'character',
	'outcome','o',1,'character',
	'outdir','1',1,'character',
	'srsid','2',1,'character',
	'sbeta','3',1,'character',
	'sstd','4',1,'character',
	'spval','5',1,'character',
	'seaf','6',1,'character',
	'threshold','t',1,'character',
	'sALLELE1','7',1,'character',
	'sALLELE0','8',1,'character',
    'ssep','9',1,'character',
    'schr','a',1,'character',
    'spos','b',1,'character',
    'sref','c',1,'character',
    'salt','d',1,'character',
	'orsid','e',1,'character',
	'obeta','f',1,'character',
	'ostd','g',1,'character',
	'opval','h',1,'character',
	'oeaf','i',1,'character',
	'oALLELE1','j',1,'character',
	'oALLELE0','k',1,'character',
    'osep','l',1,'character',
    'ochr','m',1,'character',
    'opos','n',1,'character',
    'oref','p',1,'character',
    'oalt','t',1,'character',    
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
	
Usage:
	--exposure	exposure GWAS file
	--outcome	outcome GWAS file or ieu-id
	--outdir	output dir
  	--threshold	pvalue threshold for pval
    --------------------------------------
	--rsid		rsid col id in exposure
    --chr       chr col id in exposure
    --pos       pos col id in exposure
    --ref       ref col id in exposure
    --alt       alt col id in exposure
	--beta 		beta col id in exposure
	--std		se col id in exposure
	--pval		pval col id in exposure
	--ALLELE1	allele1 col id in exposure
	--ALLELE0	allele0 col id in exposure
    --eaf       allele1 frequency in exposure
	--pop		 EUR, SAS, EAS, AFR, AMR,default [\"EUR\"]
    --sep       sep,default [\"\t\"]

    
    --help		usage
\n")
	q(status=1);
}
times<-Sys.time()
if(!is.null(opt$help)){print_usage(spec)}
if(is.null(opt$threshold)){opt$threshold = 10e-7}
if(!is.null(opt$threshold)){opt$threshold = as.numeric(opt$threshold)}
if(is.null(opt$outcome)){cat("outcome\n");print_usage(spec)}
if(is.null(opt$outdir)){opt$outdir="./"}
if(is.null(opt$srsid) & (is.null(opt$schr) | is.null(opt$spos) | is.null(opt$sref) | is.null(opt$salt))){
    cat("rsid or chr/pos/ref/alt\n")
    print_usage(spec)
}
if(is.null(opt$orsid) & (is.null(opt$ochr) | is.null(opt$opos) | is.null(opt$oref) | is.null(opt$oalt))){
    cat("rsid or chr/pos/ref/alt\n")
    print_usage(spec)
}
if(is.null(opt$std)){cat("std\n");print_usage(spec)}
if(is.null(opt$beta)){cat("beta\n");print_usage(spec)}
if(is.null(opt$pval)){cat("pval\n");print_usage(spec)}
if(is.null(opt$ALLELE1)){cat("ALLELE1\n");print_usage(spec)}
if(is.null(opt$ALLELE0)){cat("ALLELE0\n");print_usage(spec)}
if(is.null(opt$eaf)){cat("eaf\n");print_usage(spec)}
if(is.null(opt$pop)){opt$pop="EUR"}
if(is.null(opt$sep)){opt$sep='\t'}



library(TwoSampleMR)
library(ggplot2)
library(dplyr)
save.image("Rdata")
###SNP CHR BP ALLELE1 ALLELE0 A1FREQ INFO BETA_INSOMNIA SE_INSOMNIA P_INSOMNIA

###
if(is.null(opt$srsid)){
    opt$srsid="SNP"
    d<-read.table(opt$exposure,head=T)
    d=d[d[[opt$spval]] < opt$threshold,]
    db<-read.table("~/app/database/wgs_v4_download/dbsnp_156/dbsnp_156.info")
    colnames(db)=c("chr","pos","id")
    db$chr
    d$chr=as.character(d$chr)
    db$chr=as.character(db$chr)
    nd<-left_join(d,db,by=c("chr","pos"))
    nd$id[is.na(nd$id)]=paste(paste(nd$chr[is.na(nd$id)],nd$pos[is.na(nd$id)],sep=":"),nd$ref[is.na(nd$id)],nd$alt[is.na(nd$id)],sep="_")
    write.table(nd,file="tmp.e.table",sep="\t",row.names=F,quote=F)
    opt$exposure="tmp.e.table"
    opt$srsid="id"
}

escaptime=Sys.time()-times;
print(escaptime)
