library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'result','r',0,'character',
	'threshold','t',0,'character',
	'chr','c',0,'character',
	'output','o',0,'character',
	'pcol','p',0,'character',
	'lcol','l',0,'character',
	'qcol','q',0,'character',
	'xlab','x',0,'character',
	'ylab','y',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--result	 the input snp file
	--chr	the scaffold id to draw
	--threshold	the threshold default 0.9995
	--mutmap	the mutmap 
	--output	the out file 
	--pcol		the col line for point
	--lcol		the col line for line
	--qcol		the pcol line for threshold
	--xlab	
	--ylab
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$result)) { print_usage(spec)}
if ( is.null(opt$output)){ opt$output="manhattan"}
#if(is.null(opt$pcol)){opt$pcol = "delta"}
#if(is.null(opt$lcol)){opt$lcol = "slidingD"}
#if(is.null(opt$qcol)){opt$qcol = "CI"}
#if(is.null(opt$xlab)){opt$xlab = "chromosome"}
#if(is.null(opt$ylab)){opt$ylab="delta-index"}
#if(is.null(opt$thresold)){opt$threshold = 0.9995}
opt$threshold=as.numeric(opt$threshold)
times<-Sys.time()
library(ggplot2)
library(dplyr)
chr<-read.table(opt$chr);
pos<-read.table(opt$result,head=TRUE,comment.char="^")
if("pos1" %in% colnames(pos)){
	pos$pos=(pos$pos1+pos$pos2)
}
dfpos<-data.frame(CHR=pos$X.chr,BP=pos$pos,Delta=pos[[opt$pcol]]) 
if(!is.null(opt$lcol)){dfpos$Slide = pos[[opt$lcol]] }
 if(!is.null(opt$qcol)){dfpos$CI = pos[[opt$qcol]]}
dfpos = dfpos %>% filter(CHR %in% chr$V1)
#Gprime=pos$Gprime,G=pos$G,ED=pos$ED,EDprime=pos$EDprime,Gfdr=pos$Gfdr,EDfdr=pos$EDfdr)
	lev<-NULL
	lev$CHR<-levels(as.factor(dfpos$CHR))
	lev$order<-gsub("\\D","",lev$CHR)
	lev$order=as.numeric(lev$order)
	
	dfpos=merge(dfpos,lev,by="CHR")
	dfpos=arrange(dfpos,order,BP)

	dpos <- dfpos %>% group_by(order) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
	  left_join(dfpos, ., by=c("order"="order")) %>%
	  arrange(order, BP) %>%
	  mutate( BPcum=BP+tot)
	  	print("haha")
        
		axisdf <- dpos %>% group_by(CHR) %>% summarize(center=(as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )
		p1 <- ggplot(dpos) +geom_point(aes(x=BPcum, y=Delta,color=as.factor(order)))
			if(!is.null(opt$lcol)){
				p1=p1+geom_line(mapping = aes(x=BPcum,y=Slide),color="black")
			}
			if(!is.null(opt$qcol)){
				p1=p1+geom_line(mapping = aes(x=BPcum,y=CI),color="red")
			}

			if(!is.null(opt$threshold)){
				quantile= quantile(dfpos$Delta,opt$threshold,na.rm=T)
				p1=p1+geom_hline(yintercept =quantile,linetype="dashed",col="red")
			}
		    p1=p1 +
		    scale_color_manual(values = rep(c("#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2"), length(chr$V1)))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +xlab(opt$xlab)+ylab(opt$ylab)+
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )

		ggsave(file=paste(opt$out,"index.pdf",sep="."),p1,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$out,"index.png",sep="."),p1,device="png",dpi=300,height=9,width=16)
		p1 <- ggplot(dfpos) +geom_point(aes(x=BP, y=Delta,color=as.factor(order)))
		if(!is.null(opt$lcol)){
			p1=p1+geom_line(mapping = aes(x=BP,y=Slide),color="black")
		}
		if(!is.null(opt$qcol)){
			p1=p1+geom_line(mapping = aes(x=BP,y=CI),color="red")
		}
		if(!is.null(opt$threshold)){
			quantile= quantile(dfpos$Delta,opt$threshold,na.rm=T)
			p1=p1+geom_hline(yintercept =quantile,linetype="dashed",col="red")
		}
		 p1=p1+
		    theme_bw() +xlab("bp")+ylab(opt$ylab)+facet_wrap(~CHR,ncol=5)
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )

		ggsave(file=paste(opt$out,"chr.index.pdf",sep="."),p1,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$out,"chr.index.png",sep="."),p1,device="png",dpi=300,height=9,width=16)


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
