library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-m", "--merge"), type = "character", action = "store", default = NULL,help = "Input merge_result file"),
  make_option(c("-t", "--threshold"), type = "integer", action = "store",  default = 0.05 , help = "threshold default 0.05")

)
opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"))
if(is.null(opt$merge)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$threshold)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
times<-Sys.time()
library(ggplot2)
library(dplyr)
library(ggpubr)


dfs <- read.table(opt$merge, sep="\t",header=TRUE)
df=dfs %>% filter(avg_dxy > quantile(avg_dxy,1-opt$threshold)) %>% filter(avg_wc_fst  < quantile(avg_wc_fst,opt$threshold))
pop1=dfs$pop1[1]
pop2=dfs$pop2[1]
write.csv(file=paste(pop1,pop2,"dxy.select",sep="."),df,row.names=F,quote=F)
df=dfs
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

  p1 = ggplot(dpos) +geom_point(aes(x=BPcum, y=avg_dxy,color=as.factor(order)))
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
		    theme_bw() +ylab("Fst")+ xlab("Chromosome")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )

  pdf(paste(pop1, pop2,"fst_dxy.manhattan.pdf",sep="."),height=6,width=16)
  ggarrange(p1,p2,ncol=1,common.legend=F)
  dev.off()
  png(paste(pop1, pop2,"fst_dxy.manhattan.png",sep="."),height=600,width=1600)
  ggarrange(p1,p2,ncol=1,common.legend=F)
  dev.off()

  df$pi2[df$pi2==0]=min(df$pi2[df$pi2!=0])/10
  p=ggplot(dpos)+geom_point(aes(x=avg_dxy,y=avg_wc_fst,col=pi1/pi2))+xlab("Dxy")+ylab("Fst")
  ggsave(p,file=paste(pop1, pop2,"fst_dxy.scatter.png",sep="."),device="png")
  ggsave(p,file=paste(pop1, pop2,"fst_dxy.scatter.pdf",sep="."),device="pdf")
escaptime=Sys.time()-times
print("Done!")
print(escaptime)
