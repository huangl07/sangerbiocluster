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
  'window', 'r', 0, 'integer',
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

if ( is.null(opt$window))  { print_usage(spec) }
if ( is.null(opt$pvalue))  { print_usage(spec) }

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
data=output_table
colnames(data)[1:2]=c("X.chr","pos")
print(head(data))
chrlist=levels(as.factor(data$X.chr))
library(caret)

grid <- expand.grid(span = seq(0.5, 0.9, len = 10), degree = 1)
ctrl=trainControl(method="cv",number=5)

loess0<-function(chr){
   df1 =data[data$X.chr == chr,]
   if(nrow(df1) < 50){
        return=data.frame(chr=df1$X.chr,pos=df1$pos,loess=1)
   }else{
		model <- train(p ~ pos, data = df1, method = "gamLoess", tuneGrid=grid, trControl = ctrl)
        result=loess(p ~ pos,data=df1,span=model$bestTune$span) %>% predict(df1$pos)
        return=data.frame(chr=df1$X.chr,pos=df1$pos,loess=result)
   }
   return
}
rlist=lapply(chrlist,loess0)
names(rlist)=chrlist
df=do.call(rbind,rlist)
colnames(df)[1:2]=c("X.chr","pos")
df=left_join(df,output_table,by=c("X.chr","pos"))
df$logP=log10(df$loess) * -1
df$CI=opt$pvalue/nrow(df)
###输出结果2
colnames(df)[1:2]=c("X.chr","pos")
write.table(df, file = paste(opt$out,"denoise.result",sep="."), sep = "\t", quote=FALSE, row.names=FALSE)

escaptime <- Sys.time()-times
print("Done!")
print(escaptime)
q()
