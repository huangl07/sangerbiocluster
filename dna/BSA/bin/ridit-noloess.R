#!/usr/bin/env Rscript
times <- Sys.time()

if (!require("pacman")){
  install.packages("pacman")
}
library(pacman)
p_load(getopt, tidyverse, data.table, purrr)
options(scipen=200) #取消科学计算

###传参信息
spec <- matrix(c(
  'index', 'i', 0, 'character',
  'group', 's', 0, 'character',
  'out', 'o', 0, 'character',
  'pvalue', 'p', 0, 'double',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--index    the input pop.table file
	--group   the sample_id file 
	--out   the result file
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$index))   { print_usage(spec)}
if ( is.null(opt$group))  { print_usage(spec) }

if ( is.null(opt$pvalue))  { print_usage(spec) }
opt$pvalue=as.numeric(opt$pvalue)
###读取pop.table文件和sample_id文件
pop_table_info <- fread(opt$index) %>% 
  tibble()
sample_id_info <- read.table(opt$group) 
bulklist <- sample_id_info %>% filter(str_detect(V2,"B")) ## 有个bug,需要参与计算的只是bulk

print(bulklist)
###建立ridit检验数据集
select_col_names <- bulklist$V1 %>%  paste0(seq = ".AD")
df <- pop_table_info %>% select(CHROM,POS,all_of(select_col_names))
for (i in select_col_names){
    df <- df %>% separate(col=i ,into=paste(i,c("_ref","_snp"),sep=""),sep=",")
}

###ridit检验函数
my_ridit<-function(a,b){
  ref=a+b
  R=(ref/2+c(0,cumsum(ref))[1:length(ref)])/sum(ref)
  Ra=sum(R*a)/sum(a)
  Rb=sum(R*b)/sum(b)
  SR=sqrt((sum(ref*R*R)-sum(ref*R)^2/sum(ref))/(sum(ref)-1))
  z=abs(Ra-Rb)/(SR*sqrt(1/sum(a)+1/sum(b)))
  pnorm(z,lower.tail=F)
}

###ridit检验
data_ref <- select(df,paste0(select_col_names,sep='_ref')) %>% sapply(as.numeric)
data_snp <- select(df,paste0(select_col_names,sep='_snp')) %>% sapply(as.numeric)
p_value <- map(1:length(data_ref[ ,1]), function(x){my_ridit(data_ref[x, ],data_snp[x, ])}) %>% unlist()
pop_table_info$p <- p_value

###输出结果

output_table <- pop_table_info %>% 
  select(CHROM, POS, Ref, Alt, Vtype, p)



output_table=na.omit(output_table)
write.table(output_table, file = paste(opt$out,"detail.result",sep="."), sep = "\t", quote=FALSE, row.names=FALSE)

###denoise
df=output_table
colnames(df)[1:2]=c("X.chr","pos")

#ndf=df%>% group_by(X.chr,pos1,pos2) %>% summarise(meanP=mean(p),ncount=n())
ndf=df
ndf$logP=log10(ndf$p)*-1
ndf$CI=log10(opt$pvalue/nrow(ndf)) * -1 
###输出结果2
write.table(ndf, file = paste(opt$out,"denoise.result",sep="."), sep = "\t", quote=FALSE, row.names=FALSE)

escaptime <- Sys.time()-times
print("Done!")
print(escaptime)
q()
