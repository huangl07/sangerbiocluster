library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'Rdata','m',1,'character',
	'out','o',1,'character',
	'popt','t',1,'character',
    'lod','l',1,'character',
    'loc','c',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript Rqtl_CP.r --map --loc --trt --out --num

Usage:
    Rdata
    out
    popt
    lod
    loc
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggplot2)

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$Rdata) ) { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$popt) ) { print_usage(spec) }
if ( is.null(opt$lod) ) {print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec)}
if (is.null(opt$method)){opt$method="cim"}
if ( is.null(opt$pvalue) ) { opt$pvalue=0.05}
if(is.null(opt$btl)){opt$btl = F }else{opt$btl=T}
load(opt$Rdata)
markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)
if(nrow(scan.result$lod) ==0){scan.result<-summary(scan,format="tabByCol",threshold=opt$lod;pm.result=opt$lod;detail$threshold=opt$lod}

if(nrow(scan.result$lod) != 0){
	markerid1<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.high)
	markerid2<-find.marker(d,scan.result$lod$chr,scan.result$lod$ci.low)
	qtls<-makeqtl(d,chr=scan.result$lod$chr,pos=scan.result$lod$pos)
	if(length(markerid1) > 1){
	qtls<-refineqtl(cross=d,qtl=qtls,pheno.col=phe.name)
	}
	fit= try(fitqtl(cross=d,qtl=qtls,pheno.col=phe.name,get.est=TRUE),silent=TRUE)
	if("try-error" %in% class(fit)){
		fitqtl<-fitqtl(cross=d,qtl=qtls,pheno.col=phe.name)
	}else{
		fitqtl<-fitqtl(cross=d,qtl=qtls,pheno.col=phe.name,get.est=TRUE)
	}
	var=fitqtl$result.full["Model","%var"]
	if(length(qtls$name)  > 1){var=fitqtl$result.drop[,"%var"]}
	result=data.frame(data.frame(scan.result$lod),pos1=markerid1,pos2=markerid2,var=var, qname=phe.name)
	result$peak=rownames(result)
    result=left_join(result,detail,by=c(peak="marker","chr","pos","lod"))
	write.table(result,file=paste(phe.name,"qtl-result.result",sep="."),sep="\t",row.names=F,quote=F)
	csv<-read.table("../total.rqtl.csv",skip=4)
	tcsv=t(csv[csv$V1 %in% result$peak,c(1,4:ncol(csv))])
	colnames(tcsv)=tcsv[1,]
 	draw=data.frame(tcsv[-1,,drop=FALSE],trait=d$pheno[[phe.name]],check.names=FALSE)
	colnames(draw)
 	drawpxg<-function(x,drawd){
		draw=data.frame(Genotype=drawd[[x]],Pheno=drawd$trait)
		draw=na.omit(draw)
        summary=summary(factor(draw$Genotype))
		summary=summary[summary > 2]
		if(nlevels(factor(draw$Genotype)) == 1| length(summary)==1){
			return(NULL)
		}else{
			compare=combn(names(summary[summary > 2]),m=2,simplify=F)
			print(compare)
			if(opt$btl){
				p=ggplot(draw)+geom_bar(aes(x=Genotype,y=as.character(Pheno),fill=as.character(Pheno)),stat="identity")+scale_fill_npg()+labs(y="Pheno",fill = "Pheno" )
				p=p+theme_bw()
				ggsave(file=paste(phe.name,x,"pxg.pdf",sep="."),p)
            	ggsave(file=paste(phe.name,x,"pxg.png",sep="."),p)
				dtest=draw %>% group_by(Genotype,Pheno) %>% summarise(count=n())
				df=data.frame(x,p.value=fisher.test(matrix(dtest$count,nrow=2)$pvalue))			
				df$group1num=summary[df$group1]
				df$group2num=summary[df$group2]
				df
			}else{
				p=ggboxplot(draw,x="Genotype",y="Pheno",group="Genotype",fill="Genotype")+stat_compare_means(comparisons=compare,aes(label=..p.format..),method="t.test")+scale_fill_npg()
				p=p+theme_bw()
				ggsave(file=paste(phe.name,x,"pxg.pdf",sep="."),p)
            	ggsave(file=paste(phe.name,x,"pxg.png",sep="."),p)
				compare_means(data=draw[draw$Genotype %in% unlist(compare),],Pheno~Genotype,method="t.test")
				df=compare_means(data=draw[draw$Genotype %in% unlist(compare),],Pheno~Genotype,method="t.test")
				df$group1num=summary[df$group1]
				df$group2num=summary[df$group2]
				df
			}
			
		}
	}
	compare=lapply(result$peak,drawpxg,drawd=draw);
	compare_df=do.call(rbind,compare)
	write.table(file=paste(phe.name,"peaks.pxg.result",sep="."),compare_df,sep="\t",row.names=F)
	print("haha4");
}else{
	cat("",file=paste(phe.name,"qtl-result.result",sep="."))
}
	legend=pm.result
	scan=scan[order(as.numeric(as.character(scan$chr)),as.numeric(scan$pos)),]
	legend=round(legend,2)
	pdf(file=paste(phe.name,"scan.pdf",sep="."),height=9,width=16)
	plot(scan)
	print("haha")
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
	png(file=paste(phe.name,"scan.png",sep="."),height=900,width=1600)
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	legend("topright",legend=legend,col=rainbow(length(pm.result)),pch=1)
	dev.off()
	#compare=list(c("AC","BC"),c("AC","BD"),c("AD","BC"),c("AD","BD"),c("AC","AD"),)
	#drawpxg<-function(x,result,d,phe.name){
	#	draw=data.frame(Genotype=d$geno[[result$lod.chr[x]]]$data[,result$peak[x]],Pheno=d$pheno[[phe.name]])
	#	print(head(draw))
	#	draw$Genotype[draw$Genotype == "1"] = "AA"
	#	draw$Genotype[draw$Genotype == "2"] = "AB"
	#	draw$Genotype[draw$Genotype == "3"] = "BB"
	#	p=ggboxplot(draw,x="Genotype",y="Pheno",group="Genotype",fill="Genotype")+stat_compare_means(comparisons=compare,aes(label=..p.format..),method="t.test")+scale_fill_npg()
	#	p=p+theme_bw()
	#	ggsave(file=paste(result$peak[x],"pxg.pdf",sep="."),p)
	#	compare_means(data=draw,Pheno~Genotype,method="t.test")
	#}
	#compare=lapply(1:nrow(result),drawpxg,result=result,d=d,phe.name=phe.name);
	#compare_df=do.call(rbind,compare)
	#write.table("peaks.pxg.result",compare_df,sep="\t",row.names=F)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)


