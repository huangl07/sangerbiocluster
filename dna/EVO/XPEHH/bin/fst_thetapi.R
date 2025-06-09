library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
  'merge','m',1,'character',
  'threshold','t',1,'numeric',
  'help','h',0,'logical'
), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
Usage:
	--merge	merge_result file
	--thre  threshold default 0.05
	--help		usage
\n")
  q(status=1);
}
now_times <-Sys.time()
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$merge) ) { print_usage(spec) }
if ( is.null(opt$threshold)){opt$threshold=0.05}else{opt$threshold=as.numeric(opt$threshold)}

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(grid))
suppressMessages(library("gridExtra"))
suppressMessages(library(egg))
library(ggplot2)
library(egg)
df<-read.table(opt$merge,head=T)
df$pi2[df$pi2==0]=min(df$pi2[df$pi2!=0])/10
df$theta=df$pi1/df$pi2
df$fst=df$avg_wc_fst
thre1=opt$threshold
thre2=1-opt$threshold
diff<-(max(df$theta)-min(df$theta))/2000
theta1<-ceiling((quantile(df$theta,thre1,na.rm=T)-quantile(df$theta,0,na.rm=T))/diff)
print(theta1)
theta2<-ceiling((quantile(df$theta,1,na.rm=T)-quantile(df$theta,thre2,na.rm=T))/diff)
g1=ggplot(df)+theme_bw()+geom_histogram(aes(theta,y=..density..),col=c(rep("blue", theta1),rep("grey",2001-theta1-theta2),rep("green",theta2)),binwidth=diff)
g1<-g1+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+theme(panel.background=element_blank())+ylab("Theta Pi %")
g1<-g1+geom_vline(aes(xintercept=quantile(df$theta,probs=thre2,na.rm=T)),linetype=5,col="black")
g1<-g1+geom_vline(aes(xintercept=quantile(df$theta,probs=thre1,na.rm=T)),linetype=5,col="black")
empty<-ggplot()+theme(panel.background=element_blank())
scatter<-ggplot()+theme_bw()
scatter<-scatter+geom_point(aes(df$theta[df$theta > quantile(df$theta,probs=thre1) & df$theta < quantile(df$theta,probs=thre2)],df$fst[df$theta > quantile(df$theta,probs=thre1) & df$theta < quantile(df$theta,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(df$theta[df$theta < quantile(df$theta,probs=thre1) & df$fst < quantile(df$fst,probs=thre2)],df$fst[df$theta < quantile(df$theta,probs=thre1) & df$fst < quantile(df$fst,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(df$theta[df$theta < quantile(df$theta,probs=thre1) & df$fst > quantile(df$fst,probs=thre2)],df$fst[df$theta < quantile(df$theta,probs=thre1) & df$fst > quantile(df$fst,probs=thre2)]),colour='blue',fill='blue')
scatter<-scatter+geom_point(aes(df$theta[df$theta > quantile(df$theta,probs=thre2) & df$fst < quantile(df$fst,probs=thre2)],df$fst[df$theta > quantile(df$theta,probs=thre2) & df$fst < quantile(df$fst,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(df$theta[df$theta > quantile(df$theta,probs=thre2) & df$fst > quantile(df$fst,probs=thre2)],df$fst[df$theta > quantile(df$theta,probs=thre2) & df$fst > quantile(df$fst,probs=thre2)]),colour='green',fill='green')
scatter<-scatter+theme(panel.background=element_blank())+xlab("Theta Pi %")+ylab("Fst")
scatter<-scatter+geom_hline(aes(yintercept=quantile(df$fst,probs=thre2)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(df$theta,probs=thre2)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(df$theta,probs=thre1)),linetype=5,col="black")
difffst<-(max(df$fst)-min(df$fst))/2000
fst1<-ceiling((quantile(df$fst,1,na.rm=T)-quantile(df$fst,thre2,na.rm=T))/difffst)
g3<-ggplot()+theme_bw()
g3<-g3+geom_histogram(aes(df$fst,y=..density..),colour=c(rep('grey',2001-fst1),rep('orange',fst1)),binwidth=difffst)+ylim(0,8)
g3<-g3+theme(axis.title.y=element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(panel.background=element_blank())+coord_flip()+xlab("Fst distribution")
g3<-g3+geom_vline(aes(xintercept=quantile(df$fst,probs=thre2)),linetype=5,col="black")+ylab("Fst Frequency %")
  
  pdf(paste(df$pop1[1], df$pop2[1],"fst_pi.pdf",sep="."))
  ggarrange(g1,empty,scatter,g3,widths=c(3,1),heights=c(1,3))
  dev.off()
  png(paste(df$pop1[1], df$pop2[1],"fst_pi.png",sep="."))
  ggarrange(g1,empty,scatter,g3,widths=c(3,1),heights=c(1,3))
  dev.off()

df=subset(df, (df$theta < quantile(df$theta,probs=thre1) & df$fst > quantile(df$fst, probs=thre2)) |  (df$theta > quantile(df$theta,probs=thre2) & df$fst > quantile(df$fst, probs=thre2)))
df=df %>% arrange(chromosome,window_pos_1) %>% select(pop1,pop2,chromosome,window_pos_1,window_pos_2,pi1,pi2,theta,fst)
write.table(file=paste(df$pop1[1], df$pop2[1],"fst_thetapi.select",sep="."), row.names=FALSE, quote = FALSE, sep = "\t", df)


escaptime <- Sys.time()-now_times
print("Done!")
print(escaptime)
