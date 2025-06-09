rm(list = ls())
times<-Sys.time()

library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outdir','o',0,'character',
  'winsize','w',0,'character',
  'pvalue','p',0,'character',
  'LOD','l',0,'character'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--infile	the input table file
	--outdir	the output dir 
  --winsize the winsize
  --pvalue  the output pvalue threshold
  --LOD the LOD threshold
	--help		usage
\n")
	q(status=1);
}

if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outdir))  { print_usage(spec) }

if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}

library("doParallel")
library("plyr")
library("dplyr")
library("rootSolve")
library("ggplot2")
library("data.table")
if( is.null(opt$pvalue)){opt$pvalue=0.001}
if( is.null(opt$LOD)){opt$LOD=5}
if(is.null(opt$winsize)){opt$winsize=1000000}
opt$LOD=as.numeric(opt$LOD)
opt$pvalue=as.numeric(opt$pvalue)
LOD<-function(dir= NULL,filegen = NULL,width = NULL){
  df<-read.csv(filegen,sep="\t",header = TRUE,stringsAsFactors=FALSE)

  data1=data.frame(df$X.chr,pos=df$pos,A1=df$n4,A2=df$n2,a1=df$n3,a2=df$n1)
  Resolution="HIGH"
  Plotformat1="png"
  Calculate<-function(data,width){
    # columns are:
    # 1. chromosome
    # 2. chromosome coordinate 
    # 3. observed number of allele A in the pool with low phenotype
    # 4. observed number of allele a in the pool with low phenotype
    # 5. observed number of allele A in the pool with high phenotype
    # 6. observed number of allele a in the pool with high phenotype
    chr = data1[,1];cc<-unique(chr)
    ff<-numeric();tt<-numeric()
    for(iii in 1:length(cc)){
      data<-data1[which(data1[,1]==cc[iii]),]
      if(nrow(data) < 50){next;}
      chrom = data[,1]; pos=data[,2];e = data[,3];f = data[,4];g = data[,5];h = data[,6]
      lfrq<-cbind(e,f);hfrq<-cbind(g,h)
      Al<-lfrq[,1];al<-lfrq[,2];Ah<-hfrq[,1];ah<-hfrq[,2];pos<-as.matrix(pos)
      P_QTGL<-Al/(Al+al);P_QTGH<-Ah/(Ah+ah)
      P_QTGL<-matrix(P_QTGL,ncol = 1);P_QTGH<-matrix(P_QTGH,ncol = 1)
      nn<-dim(P_QTGH)[1]
      nQTG<-cbind(Al,Ah)
      lod<-matrix(NA,nrow = nn,ncol = 1)
      for(i in 1:nn){
        c1<-Al+al;c2<-Ah+ah
        unexist_QTN<-(choose(c1[i],nQTG[i,1])*((1/2)^c1[i]))*(choose(c2[i],nQTG[i,2])*((1/2)^c2[i]))
        exist_QTN<-choose(c1[i],nQTG[i,1])*((P_QTGL[i])^nQTG[i,1])*((1-P_QTGL[i])^(c1[i]-nQTG[i,1]))*choose(c2[i],nQTG[i,2])*((P_QTGH[i])^nQTG[i,2])*((1-P_QTGH[i])^(c2[i]-nQTG[i,2]))
        lod[i]<-log10(exist_QTN/unexist_QTN)
      }
      result<-cbind(pos,lod)
      x=pos
      y=lod     
      w=width
      olen = length(x)
      #beginning intervals
      xfirst = x[1,]      #start site
      startband = x <= xfirst + w  
      xstart = xfirst - (x[startband] - xfirst)
      ystart = y[startband]
      nstart = length(xstart)-1
      #end of intervals
      xlast=x[length(x),]#last site
      endband = x >= xlast - w
      xend = xlast + (xlast - x[endband])
      yend = y[endband]
      nend = length(xend) - 1   
      a=rev(xstart)   #reversal
      a=a[-length(a)] #remove the last site of b
      b=rev(xend)     
      b=b[-1]         #remove the first site of b
      x<-c(a,x,b)
      c=rev(ystart)
      c=c[-length(c)] #remove the last site of b
      d=rev(yend)     
      d=d[-1]
      y<-c(c,y,d)
      lx = length(x);ly= length(y)

      cl<- makeCluster(6)
      registerDoParallel(cl)
      yw=foreach(i=1:lx)%dopar%
      {
        c = x[i]
        inband = ( x >= (c-w)&x <= (c+w))
        xfrac = (abs(x[inband] - c))/w
        xwt = (1.0 - abs(xfrac)^3)^3    #Weights 
        xwt[abs(xfrac) >= 1] = 0    
        ywin = sum(y[inband]*xwt)/sum(xwt)
      }
      stopCluster(cl)
      g2=cbind(y,yw)
      h1=length(a)+1
      h2=length(a)+olen
      g2=g2[h1:h2,]
      smooth_G=cbind(pos,g2)
      colnames(smooth_G)<-NULL
      G1<-as.matrix(smooth_G)
      G<-matrix(as.numeric(G1),nrow = nrow(G1))
      smoothG<-G[,3]*2*log(10)  
      smoothG<-as.matrix(smoothG,ncol=1)
      prob<-function(x){1-pchisq(x, 1)}
      p<-apply(smoothG,1,prob)
      p<-as.matrix(p,ncol=1)
      G<-cbind(chrom,G,p)
      colnames(G)<-c("Chrom","Pos","LOD","Smooth_LOD","P_Value")
      ff<-rbind(ff,G)
    }
    output=as.data.frame(ff)
    colnames(output)=c("Chrom","Pos","LOD","Smooth_LOD","P_Value")
    return(output)
  }
  calculateresult<-Calculate(data,width)
  write.table(calculateresult,paste(dir,"/","tmp.csv",sep=""),sep="\t",row.names=F,col.names = T,quote=F) 
  #write.table(calculateresult$result,paste(dir,"/","result.csv",sep=""),sep=",",row.names=F,quote=F,col.names = T) 
}

draw<-function(drawdata=NULL,col=NULL){
      x=data.frame(drawdata)
      colnames(x)<-c("chr","bp","lod","smoothlod","logp")
      #print(head(x))
      x$CHR=as.numeric(as.factor(x$chr))
      CHR=BP=P=index=NULL
      d=data.frame(CHR=x$CHR, BP=x$bp, lod=x$lod,smoothlod=x$smoothlod,logp=x$logp)
      d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(lod) & is.numeric(smoothlod)))
      d <- d[order(d$CHR, d$BP), ]
      d$pos=NA
      d$index=NA
      ind = 0
      for (i in unique(d$CHR)){
        ind = ind + 1
        d[d$CHR==i,]$index = ind
      }

      nchr = length(unique(d$CHR))
      if (nchr==1) { ## For a single chromosome
        d$pos=d$BP/1000000
        xlabel = paste('Chromosome',unique(d$CHR),'(Mb)')
      } else { ## For multiple chromosomes
        lastbase=0
        ticks=NULL
        for (i in unique(d$index)) {
          if (i==1) {
            d[d$index==i, ]$pos=d[d$index==i, ]$BP
          } else {
            lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
            d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
          }
          ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
        }
        xlabel = 'Chromosome'
        labs <- unique(d$CHR)
      }
      
      
      # Initialize plot
      xmax = ceiling(max(d$pos) * 1.03)
      xmin = floor(max(d$pos) * -0.03)
      ymax<-(ceiling(max(d$lod)*1.25/10))*10
      margin_space<-0.5
      par(mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
      suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                                        xlim=c(xmin,xmax),ylim=c(0,ymax),
                                        xlab="", cex.lab=0.6,ylab=""))
      dotargs <- as.list(match.call())[-1L]
      dotargs <- list()
      ## And call the plot function passing NA, your arguments, and the default
      do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
      # Create a vector of alternatiting colors
      col=rep(col, max(d$CHR))
      
      
      
      # Add points to the plot
      if (nchr==1) {
        with(d,  points(pos, lod, col="LightGray", pch=20,cex=0.2))
      } else {
        # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
        icol=1
        for (i in unique(d$index)) {
          #with(d[d$index==unique(d$index)[i], ], points(pos, lod, col="LightGray", pch=20,cex=0.1))
          with(d[d$index==unique(d$index)[i], ], lines(pos, smoothlod, col="darkGray", pch=20,cex=0.1))
          icol=icol+1
        }
      }
      suppressWarnings(axis(4,ylim=c(0,ymax),cex.lab=0.5,mgp=c(3,-0.2,0),tcl=-0.2,cex.axis=0.4,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray"))
      abline(h=opt$LOD,col="red",lty=2)
      mtext("smoothLOD",side=4,line=0.5,cex=0.5,font = 1,col="DarkGray")
      
      
      
      par(new=TRUE,mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
      suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=20, yaxt="n" ,
                                        xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))),
                                        xlab=xlabel, ylab="",cex.lab=0.55,font = 1))#line=0.78,
      dotargs <- list()
      # And call the plot function passing NA, your arguments, and the default
      do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
      # Add an axis. 
      if (nchr==1) { #If single chromosome, ticks and labels automatic.
        suppressWarnings(axis(1,cex.lab=0.5,mgp=c(3,-0.15,0.3),tcl=-0.2,cex.axis=0.4))
      } else { # if multiple chrs, use the ticks and labels you created above.
        suppressWarnings(axis(1, at=ticks, labels=labs,cex.lab=0.5,mgp=c(3,-0.15,0.3),tcl=-0.2,cex.axis=0.4))
      }
      # Create a vector of alternatiting colors
      col=rep(col, max(d$CHR))
      # Add points to the plot
      if (nchr==1) {
        with(d, lines(pos, logp, pch=20, cex=0.1,col=col[1]))
      } else {
        # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
        icol=1
        for (i in unique(d$index)) {
          with(d[d$index==unique(d$index)[i], ], points(pos, logp, col="blue",pch=20,cex=0.18))
          icol=icol+1
        }
      }
      axis(2,ylim=c(0,max(d$logp)),cex.lab=0.5,mgp=c(3,0.2,0),tcl=-0.2,cex.axis=0.4,col.axis = col[1],col.ticks=col[1],col =col[1] )
      abline(h=log10(opt$pvalue/nrow(df))*-1,col="blue",lty=3)
      mtext("-log10(p)",side=2,line=0.8,cex=0.5,font = 1,col = col[1])
}

LOD(dir=opt$outdir,filegen = opt$infile,width = as.numeric(opt$winsize))
df<-read.table("tmp.csv",head=T)
colnames(df)=c("X.chr","pos","LOD","smooth_LOD","pvalue")
df$pvalue[df$pvalue==0]=min(df$pvalue[df$pvalue !=0])*0.1
df$logp=log10(df$pvalue)*-1
write.table(df,paste(opt$outdir,"/","result.txt",sep=""),sep="\t",row.names=F,col.names = T,quote=F) 
opt$pvalue=opt$pvalue/nrow(df)
df$SlidingD=0.9
df$SlidingD[which(df$pvalue<=opt$pvalue&df$smooth_LOD>=opt$LOD)]=1
df$CI=0.95
write.table(df,paste(opt$outdir,"/","significant.result.txt",sep=""),sep="\t",row.names=F,col.names = T,quote=F) 

Resolution="HIGH"
Plotformat1="png"
manwidth=15000; manhei=8000;units= "px";manwordre =30;manfigurere=1000
col=c("blue","red")
png(paste(opt$outdir,"/","LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
draw(df[,c(1:4,6)],col)
dev.off()
pdf(paste(opt$outdir,"/","LOD.pdf",sep=""),width=15, height=8)
draw(df[,c(1:4,6)],col)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
