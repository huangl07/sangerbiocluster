library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'exposure','e',1,'character',
	'outcome','o',1,'character',
	'outdir','d',1,'character',
	'rsid','r',1,'character',
	'beta','b',1,'character',
	'std','s',1,'character',
	'pval','p',1,'character',
	'eaf','f',1,'character',
	'threshold','t',1,'character',
	'ALLELE1','4',1,'character',
	'ALLELE0','0',1,'character',
    'sep','q',1,'character',
    'chr','c',1,'character',
    'pos','1',1,'character',
    'ref','2',1,'character',
    'alt','3',1,'character',
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
    --------------------------------------
	--rsid		rsid col id in exposure
    --------------------------------------
    --chr       chr col id in exposure
    --pos       pos col id in exposure
    --ref       ref col id in exposure
    --alt       alt col id in exposure
    --------------------------------------
	--beta 		beta col id in exposure
	--std		se col id in exposure
	--pval		pval col id in exposure
	--threshold	pvalue threshold for pval
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
if(is.null(opt$rsid) & (is.null(opt$chr) | is.null(opt$pos) | is.null(opt$ref) | is.null(opt$alt))){
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
if(is.null(opt$rsid)){
    opt$rsid="SNP"
    d<-read.table("meta_fat.txt",head=T)
    d=d[d[[opt$pval]] < opt$threshold,]
    db<-read.table("~/app/database/wgs_v4_download/dbsnp_156/dbsnp_156.info")
    colnames(db)=c("chr","pos","id")
    db$chr
    d$chr=as.character(d$chr)
    db$chr=as.character(db$chr)
    nd<-left_join(d,db,by=c("chr","pos"))
    nd$id[is.na(nd$id)]=paste(paste(nd$chr[is.na(nd$id)],nd$pos[is.na(nd$id)],sep=":"),nd$ref[is.na(nd$id)],nd$alt[is.na(nd$id)],sep="_")
    write.table(nd,file="tmp.table",sep="\t",row.names=F,quote=F)
    opt$exposure="tmp.table"
    opt$rsid="id"
}
print(opt)
exposure_data <- read_exposure_data(
  filename = opt$exposure,
  sep = opt$sep,
  snp_col = opt$rsid,
  beta_col = opt$beta,
  se_col = opt$std,
  effect_allele_col = opt$ALLELE1,
  other_allele_col = opt$ALLELE0,
  pval_col = opt$pval,
  eaf_col= opt$eaf,
)

opt$threshold=10e-5
exposure_data_clump <- clump_data(exposure_data[exposure_data$pval.exposure < opt$threshold,],clump_r2=0.01,pop = opt$pop)

if(!file.exists(opt$outcome)){
    outcome_data <- extract_outcome_data(
    snps = exposure_data_clump$SNP,
    outcomes = opt$outcome
    )
}else{
    outcome_data <- read_outcome_data(
        filename = opt$outcome,
        sep = opt$sep,
        snp_col = opt$rsid,
        beta_col = opt$beta,
        se_col = opt$std,
        effect_allele_col = opt$ALLELE1,
        other_allele_col = opt$ALLELE0,
        pval_col = opt$pval,
        eaf_col= opt$eaf,
    )
    outcome_data_clump <- clump_data(outcome_data[outcome_data$pval.outcome < opt$threshold,],clump_r2=0.01,pop = opt$pop)
    dat <- harmonise_data(
        exposure_dat = exposure_data_clump, 
        outcome_dat = outcome_data_clump
    )
}
if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}
setwd(opt$outdir)
write.table(file="mr.data",sep="\t",dat,row.names=F,quote=F)
res <- mr(dat)
hetero=mr_heterogeneity(dat)
pelio=mr_pleiotropy_test(dat)
write.table(file="mr.result",res,sep="\t",quote=F,row.names=F)
write.table(file="mr.hetero",hetero,sep="\t",quote=F,row.names=F)
write.table(file="mr.pelio",pelio,sep="\t",quote=F,row.names=F)
res_single <- mr_singlesnp(dat)
write.table(file="mr.singlesnp.result",res_single,sep="\t",quote=F,row.names=F)
p1 <-mr_scatter_plot(res, dat)
ggsave(file="mr.pdf",p1[[1]])
if(nrow(res_single))
p2 <- mr_forest_plot(res_single)
ggsave(file="mr.forest.pdf",p2[[1]])
if(nrow(res_single) > 50){
    top50=res_single[order(res_single$p,decreasing=F)[1:50],]
    p4 <- mr_forest_plot(top50)
    ggsave(file="mr.forest.top50.pdf",p4[[1]])
}
p3 <- mr_funnel_plot(res_single)
ggsave(file="mr.funnel.pdf",p3[[1]])


escaptime=Sys.time()-times;
print(escaptime)
