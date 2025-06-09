
#!/usr/bin/env Rscript
times<-Sys.time()

if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)

###传参信息
spec <- matrix(c(
  #'degfile', 'd', 0, 'character',
  'maf1', '1', 0, 'character',
  'maf2', '2', 0, 'character',
  'out', 'o', 0, 'character',
  'pvalue','p',0,'character',
  #'db','b','0','character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--maf1   输入分组的maf文件
    --maf2   输入分组的maf文件	
    --out    输出文件名字
    --pvalue p值
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$maf1))   { print_usage(spec)}
if ( is.null(opt$maf2))  { print_usage(spec) }
if ( is.null(opt$out))  { print_usage(spec) }
if(is.null(opt$pvalue)){opt$pvalue = 0.05}
library(dplyr)
LC<-read.table(opt$maf1,head=T,sep="\t",comment.char="^",quote="")
IDC<-read.table(opt$maf2,head=T,sep="\t",comment.char="^",quote="")
LCnSample=nlevels(factor(LC$Tumor_Sample_Barcode))
IDCnSample=nlevels(factor(IDC$Tumor_Sample_Barcode))
d1=LC %>% group_by(Chromosome,vcf_pos,HGVSc,Hugo_Symbol) %>% summarise(count=n(),all= LCnSample)
d2=IDC %>% group_by(Chromosome,vcf_pos,HGVSc,Hugo_Symbol) %>% summarise(count=n(),all= IDCnSample)
merged=full_join(d1,d2,by=c("Chromosome","vcf_pos","HGVSc","Hugo_Symbol")) 
merged$all.x[is.na(merged$all.x)]=LCnSample
merged$all.y[is.na(merged$all.y)]=IDCnSample
merged$count.x[is.na(merged$count.x)]=0
merged$count.y[is.na(merged$count.y)]=0
myfun=function(d,data){fisher.test(matrix(c(data$count.x[d],data$all.x[d],data$count.y[d],data$all.y[d]),nrow=2))$p.value}
merged$pvalue=unlist(lapply(1:nrow(merged),myfun,data=merged))
colname(merged)=c("Chromosome","Position","HGVSc","Hugo_Symbol","group1_allele","group1_all","group2_allele","group2_all","pvalue")
result=merged %>% filter(pvalue < as.numeric(opt$pvalue))
write.table(file=paste(opt$out,"table",sep="."),row.names=F,sep="\t",merged)
write.table(file=paste(opt$out,"signals",sep="."),row.names=F,sep="\t")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()

