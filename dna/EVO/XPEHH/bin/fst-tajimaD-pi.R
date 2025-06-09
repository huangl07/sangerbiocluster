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
	--threshold  threshold default 0.05
	--help		usage
\n")
  q(status=1);
}
now_times <-Sys.time()
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$merge) ) { print_usage(spec) }
if ( is.null(opt$threshold)){opt$threshold=0.05}else{opt$threshold=as.numeric(opt$threshold)}


library(ggplot2)
library(ggpubr)
library(dplyr)

df<-read.table(opt$merge,head=T)
df$BP=(df$window_pos_1+df$window_pos_2)/2
df$CHR=df$chromosome
lev=NULL
lev$CHR<-levels(as.factor(df$CHR))
lev$order<-gsub("\\D","",lev$CHR)
lev$order=as.numeric(lev$order)
dfpos=merge(df,lev,by="CHR")
dfpos=arrange(dfpos,order,BP)
dpos <- dfpos %>% group_by(order) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
  left_join(dfpos, ., by=c("order"="order")) %>%
  arrange(order, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- dpos %>% group_by(CHR) %>% summarize(center=(as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )
  p1 = ggplot(dpos) +geom_line(aes(x=BPcum, y=pi1),color="blue")+geom_point(aes(x=BPcum, y=pi2),color="green")
  p1=p1 + scale_color_manual(values = rep(c("#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2"), nlevels(as.factor(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +ylab("Dxy")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank(),
              axis.title.x=element_blank()
		    )
  p2 = ggplot(dpos) +geom_point(aes(x=BPcum, y=avg_wc_fst,color=as.factor(order)))
  p2=p2 + scale_color_manual(values = rep(c("#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2"), nlevels(as.factor(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +ylab("Fst")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
  p3 = ggplot(dpos) +geom_line(aes(x=BPcum, y=tajimaD1),color="blue")+geom_point(aes(x=BPcum, y=tajimaD2),color="green")
  p3=p3 + scale_color_manual(values = rep(c("#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2"), nlevels(as.factor(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +ylab("Tajima'D")+ xlab("Chromosome")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
  pdf(paste(df$pop1[1], df$pop2[1],"fst_pi_tajimaD.manhattan.pdf",sep="."),height=6,width=16)
  ggarrange(p1,p2,ncol=1,common.legend=F)
  dev.off()
  png(paste(df$pop1[1], df$pop2[1],"fst_pi_tajimaD.manhattan.png",sep="."),height=600,width=1600)
  ggarrange(p1,p2,ncol=1,common.legend=F)
  dev.off()

escaptime <- Sys.time()-now_times
print("Done!")
print(escaptime)
