#!/usr/bin/env Rscript
start_time <- Sys.time()
library(getopt)
spec <- matrix(c(
  'csv', 'c', 0, 'character',
  'out', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--csv    输入csv表格
  --out    输出基因突变全景图
	--help   usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$csv))   { print_usage(spec) }
if ( is.null(opt$out))   { print_usage(spec) }
library(ggplot2)
library(patchwork)
data <- read.csv(opt$csv,header=TRUE,sep=",")
data_postive<- data[c(1,2,11)]
names(data_postive)<-c("Gene1","Gene2","Pvalue")
data_negative<- data[c(2,1,12)]
names(data_negative)<-c("Gene1","Gene2","Pvalue")
data_negative$Pvalue=data_negative$Pvalue*(-1)
alldata<-rbind(data_postive,data_negative)
display_data <- alldata
display_data$Pvalue<-abs(display_data$Pvalue)
display_data$Pvalue[which((display_data$Pvalue<"0.05")&(display_data$Pvalue>"0"))]<-"*"
display_data$Pvalue[which((display_data$Pvalue<"0.1")&(display_data$Pvalue>="0.05"))]<-"·"
display_data$Pvalue[which(display_data$Pvalue>="0.1")]<-""
display_data$Pvalue[which(display_data$Pvalue=="0")]<-""

names(display_data)<-c("Gene1","Gene2","log10Pvalue")
data_postive$Pvalue=abs(log10(data_postive$Pvalue))
data_negative$Pvalue=data_negative$Pvalue*(-1)
data_negative$Pvalue=abs(log10(data_negative$Pvalue))*(-1)
rundata<-rbind(data_postive,data_negative)
rundata<-cbind(rundata,display_data$log10Pvalue)
names(rundata)<-c("Gene1","Gene2","Pvalue","log10Pvalue")

legend1<-data.frame(y=c(-3,-2,-1,0,1,2,3),x=c(""))

p1<-ggplot(rundata,aes(x=Gene1,y=Gene2,fill=Pvalue))+
    geom_tile(color="grey80")+
    geom_text(aes(x=Gene1, y=Gene2, label=log10Pvalue)) +
    theme(legend.position = "none", aspect.ratio=1,axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90)) +
    scale_fill_gradient2(low="#003366", high="#990033", mid="white",breaks=c(-3,-2,-1,0,1,2,3))
       
p3<-ggplot(legend1[4:7,],(aes(x=x,y=y))) +
  geom_tile(aes(fill=y),color="grey80")+
  scale_fill_gradient(low="white", high="#990033",breaks=c(0,1,2,3)) +
  theme(axis.ticks.x = element_blank(),plot.title=element_text(size=11),axis.title.y=element_text(size=11),axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(title="Co_occurance",y = "-log10(Pvalue)")+
  coord_fixed(ratio=1) +
  scale_y_continuous(position = "left",breaks=c(0,1,2,3),labels=c("0","1","2","3"))

p4<-ggplot(legend1[4:7,],(aes(x=x,y=y))) +
  geom_tile(aes(fill=y),color="grey80")+
  scale_fill_gradient(low="white", high="#003366",breaks=c(0,1,2,3)) +
  theme(axis.ticks.x = element_blank(),plot.title=element_text(size=11),axis.title.y=element_text(size=11),axis.title.x=element_text(size=11,hjust=0.1,vjust=1)) +
  guides(fill=FALSE) +
  labs(title="Exclusive",y = "-log10(Pvalue)",x="* p<0.05\n· p<0.1")+
  coord_fixed(ratio=1) +
  scale_y_continuous(position = "left",breaks=c(0,1,2,3),labels=c("0","1","2","3"))


patch<-p1+(plot_spacer()/p3/p4/plot_spacer()+plot_layout(heights = c(1,3,3,2)))+plot_layout(widths = c(5, 1))
ggsave(paste0(opt$out,"mutation_interaction_map.png"),patch,units="in",dpi=300, width=4, height=4, device="png")
ggsave(paste0(opt$out,"mutation_interaction_map.pdf"),patch,units="in", width=8, height=8, device="pdf")