library('getopt')
options(bitmapType='cairo')

spec = matrix(c(
				'help',   'h', 0, 'logical',  'print help document, and exit',
				'mark', 'm', 1, 'character','file, required',
				'outdir', 'o', 1, 'character','outdir, optional',
				'name', 'n', 1, 'character','out file prefix, optional'

		), byrow=TRUE, ncol=5);
opt=getopt(spec)
print_usage <- function(spec=NULL){
	script <- get_Rscript_filename()
	cat(getopt(spec, usage=TRUE));
	cat("Usage example :\n")
	cat("
Usage example:
        Rscript Haplotype.r -m markfile -o outdir -n name
Usage:
  -m     map   file
  -o     out  dir
  -n     prefix name
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) )	{ print_usage(spec) }
if ( is.null(opt$outdir) ) { opt$outdir=getwd()}
if (is.null(opt$name)) { opt$name="haplotype" }


library(ggplot2)
library(reshape2)
library(tidyverse)
library(forcats)
cols=c("red","blue","orange","grey40")
names(cols)=c("A","B","H","-")
df<-read.csv(opt$mark,header=T)
colnames(df)[2:3]=c("chr","pos")
dd<-melt(df,id.vars=c("Genotype","chr","pos"))%>%
  select(-Genotype)%>%group_by(chr,pos,variable)%>%
  mutate(count_n=fct_lump(value,n=1))%>%
  group_by(chr,pos,variable,value,count_n)%>%
  summarise(n())%>%
  filter(count_n!="Other")
chr_length<-group_by(df,chr)%>%
  summarise(chrlength=max(pos))%>%
  mutate(cunsum_length=cumsum(chrlength))
chr_length$chrlength1=c(0,chr_length$cunsum_length[1:nrow(chr_length)-1])
ff<-left_join(dd,chr_length,by=c("chr"="chr"))
ff$finally_length=ff$pos+ff$chrlength1
ff$length=c(1:nrow(ff))
dm=ff %>% group_by(chr)%>% summarise(at=mean(length))
p<-ggplot(data=ff,aes(x=length,y=variable,color=value))+
geom_tile()+
scale_color_manual(values=cols)+
labs(title="Haplotype",x="",y="")+
scale_x_continuous(breaks=dm$at,labels=dm$chr)+
theme_bw()+
theme(legend.position="None",
      plot.title=element_text(hjust=0.5),
      axis.text.y=element_blank(),
      axis.ticks=element_blank())
png(filename=paste(opt$outdir,"/",opt$name,".png",sep=""),height=600,width=1000)
print(p)
dev.off()

pdf(file=paste(opt$outdir,"/",opt$name,".pdf",sep=""), height = 600, width = 1000)
print(p)
dev.off()
