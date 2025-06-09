library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-x", "--xpehh"), type = "character", action = "store", help = "Input xpclr file"),
  make_option(c("-t", "--threshold"), type = "numeric", action = "store", default = 0.05, help = "threshold for select")
)
opt = parse_args(OptionParser(option_list = option_list,
                              usage = "Usage: %prog [options] \nDescription: This Script is used to draw xpclr manhattan!"))
if(is.null(opt$xpehh)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}

times<-Sys.time()

library(ggplot2)
library(dplyr)
library(ggpubr)


dfs <- read.table(opt$xpehh, sep="\t",header=TRUE)
dfs=dfs[grep("chr",dfs$id),]
colnames(dfs)[4]="normxpehh"
df=dfs %>% filter(normxpehh >quantile(as.numeric(dfs$normxpehh),1-opt$threshold,na.rm=T) |  normxpehh <quantile(as.numeric(dfs$normxpehh),1-opt$threshold,na.rm=T))
pop=gsub(".norm.result","",opt$xpehh)
print(pop)
print(nrow(df))
write.csv(file=paste(pop,"xpehh.select",sep="."),df,row.names=F,quote=F)
df=dfs
print(nrow(df))
df$normxpehh=as.numeric(df$normxpehh)
df$BP=(as.numeric(df$pos1)+as.numeric(df$pos2))/2
quantile1=quantile(df$normxpehh,1-opt$threshold,na.rm=T)
quantile2=quantile(df$normxpehh,opt$threshold,na.rm=T)
print("haha")
print(quantile1)
print(quantile2)
df$CHR=df$id
lev=NULL
lev$CHR<-levels(as.factor(df$CHR))
lev$order<-gsub("\\D","",lev$CHR)
lev$order=as.numeric(lev$order)
dfpos=merge(df,lev,by="CHR")
dfpos=arrange(dfpos,order,BP)
print("haha")
dpos <- dfpos %>% group_by(order) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
  left_join(dfpos, ., by=c("order"="order")) %>%
  arrange(order, BP) %>%
  mutate( BPcum=BP+tot)
axisdf <- dpos %>% group_by(CHR) %>% summarize(center=(as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )
dpos$col=1
dpos$col[dpos$order %%2 ==0]=2
dpos$col[dpos$normxpehh > quantile1]=3
  p1 = ggplot(dpos) +geom_point(aes(x=BPcum, y=normxpehh,color=as.factor(col)))+geom_hline(yintercept =quantile1,linetype="dashed",col="red")
  p1=p1 + scale_color_manual(values = rep(c("slategray3", "skyblue4","red"), nlevels(as.factor(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +ylab("XP-EHH")+
		    theme( 
		      legend.position="none",
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank(),
          axis.text.x=element_text(angle = 90 , vjust = 0.5 , hjust= 1 ),
		    )
  pop=gsub(".xpehh.result","",opt$xpehh)
  ggsave(p1,file=paste(pop,"xpehh.manhattan.png",sep="."),device="png",height=9,width=16)
  ggsave(p1,file=paste(pop,"xpehh.manhattan.pdf",sep="."),device="pdf",height=9,width=16)

escaptime=Sys.time()-times
print("Done!")
print(escaptime)

escaptime=Sys.time()-times
print("Done!")
print(escaptime)
