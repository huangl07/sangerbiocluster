library('getopt');

options(bitmapType='cairo')
spec = matrix(c(
    'map','m',1,'character',
    'loc','l',1,'character',
    'trt','t',1,'character',
    'out','o',1,'character',
    'num','n',1,'character',
    'lod','d',1,'character',
    'btl','b',1,'character',
    'pop','p',1,'character',
    'method','x',1,'character',
    'pvalue','v',1,'character',
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
    --map    map file
    --loc    loc file
    --trt    trt file
    --out    out dir
    --num    pm number
    --btl
    --pvalue    select pvalue
    --lod    threshold lod value
    --method default cim
    --help        usage
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
if ( is.null(opt$map) ) { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }else{opt$num=as.numeric(opt$num)}
if ( is.null(opt$out) ) { opt$out="./";}
if ( is.null(opt$lod) && is.null(opt$pvalue) ) {print("ERROR qtl paramter");quit(status=1);}

if(is.null(opt$method)){opt$method="cim"}
if(!dir.exists(opt$out)){dir.create(opt$out)}
if(is.null(opt$btl)){opt$btl = F }else{opt$btl=T}

if(grepl("RIL",opt$pop) ){opt$pop="riself"}
d<-read.cross(file=opt$loc,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="*",genotypes=c("a","b","H"))
csv<-read.csv(opt$loc,head=F)
csv=csv[-2][-2]
if(!dir.exists(opt$out)){dir.create(opt$out)}
setwd(opt$out);
set.seed(12345)
d<-jittermap(d)
d<-calc.genoprob(d)
d=sim.geno(d)
phe.name<-colnames(d$pheno)[2]
d$pheno[[phe.name]]=as.numeric(d$pheno[[phe.name]])
pdf(paste(phe.name,"pheno.pdf",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
png(paste(phe.name,"pheno.png",sep="."))
plotPheno(d,pheno.col=phe.name)
dev.off()
model="normal"
effectscan=effectscan(d,pheno.col=phe.name)
if(opt$btl){model="binary";opt$method="scanone"}
if(opt$method == "cim"){
    scan<-cim(d,pheno.col=phe.name);
    if(is.null(opt$lod)){
        scan.pm<-cim(d,pheno.col=phe.name,n.perm=opt$num);
    }
    if(median(summary(scan)$lod) > 20){
        scan<-scanone(d,pheno.col=phe.name,model=model);
        if(is.null(opt$lod)){
            scan.pm<-scanone(d,pheno.col=phe.name,n.perm=opt$num,model=model);
        }
    }
}else{
    scan<-scanone(d,pheno.col=phe.name,model=model);
    if(is.null(opt$lod)){
        scan.pm<-scanone(d,pheno.col=phe.name,n.perm=opt$num,model=model);
    }
}
if(!is.null(opt$pvalue)){
    opt$pvalue=as.numeric(opt$pvalue)
    qtl.result=summary(scan,format="tabByCol",perms=scan.pm,alpha=opt$pvalue,drop=1)
    pm.result=summary(scan.pm,alpha=opt$pvalue)
}else{
    qtl.result<-summary(scan,format="tabByCol",threshold=as.numeric(opt$lod),drop=1)
    pm.result<-as.numeric(opt$lod)
}
markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)
detail=data.frame(marker=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod,pm=as.numeric(pm.result))
effectscan$marker=rownames(effectscan)
detail=left_join(detail,effectscan,by=c("marker","chr","pos"))
print(qtl.result)
save.image("../qtl.Rdata")
if(nrow(qtl.result$lod) ==0){qtl.result<-summary(scan,format="tabByCol",threshold=max(detail$lod)-0.5,drop=1);pm.result=max(detail$lod)-0.5;detail$threshold=max(detail$lod)-0.5}
write.table(detail,file=paste(phe.name,"detail.result",sep="."),sep="\t",row.names=F,quote=F)

if(nrow(qtl.result$lod)!=0){
    markerid1<-find.marker(d,qtl.result$lod$chr,qtl.result$lod$ci.high)
    markerid2<-find.marker(d,qtl.result$lod$chr,qtl.result$lod$ci.low)
    ids=find.marker(d,chr=qtl.result$lod$chr,pos=qtl.result$lod$pos)
    qtls<-makeqtl(d,chr=qtl.result$lod$chr,pos=qtl.result$lod$pos)
    if(length(markerid1) > 1){
    qtl<-refineqtl(cross=d,qtl=qtls,pheno.col=phe.name)
    }
    png("effect.png",units="in",width=12,height=8,res=300)
    fit= try(fitqtl(cross=d,qtl=qtls,pheno.col=phe.name,get.est=TRUE),silent=TRUE)
    if("try-error" %in% class(fit)){
        fitqtl<-fitqtl(cross=d,qtl=qtls,pheno.col=phe.name)
    }else{
        fitqtl<-fitqtl(cross=d,qtl=qtls,pheno.col=phe.name,get.est=TRUE)
    }
    dev.off()
    print("haha")
    if(length(markerid1) == 1){
        var=fitqtl$result.full["Model","%var"]
    }else{

        var=fitqtl$result.drop[,"%var"]
    }
    result=data.frame(data.frame(qtl.result$lod),pos1=markerid1,pos2=markerid2,var=var,qname=phe.name)
    nqtl=nrow(result)
    result$peak=rownames(result)
    result=left_join(result,detail,by=c(peak="marker","chr","pos","lod"))
    write.table(result,file=paste(phe.name,"qtl-result.result",sep="."),sep="\t",row.names=F,quote=F)
    if(is.null(opt$lod)){
        write.table(file=paste(phe.name,".pm.csv",sep=""),sep="\t",scan.pm,quote=F,row.names=T);
    }
    drawpxg<-function(x,result,d,phe.name){
        draw=data.frame(Genotype=d$geno[[result$chr[x]]]$data[,result$peak[x]],Pheno=d$pheno[[phe.name]])
        draw=na.omit(draw)
        print(head(draw))
        if(opt$pop != "dh"){
            draw$Genotype[draw$Genotype == "1"] = "AA"
            draw$Genotype[draw$Genotype == "2"] = "AB"
            draw$Genotype[draw$Genotype == "3"] = "BB"
        }else{
            draw$Genotype[draw$Genotype == "1"] = "AA"
            draw$Genotype[draw$Genotype == "2"] = "BB"
        }
        draw$Genotype=factor(draw$Genotype, levels=c("AA","AB","BB"))
        summary=summary(factor(draw$Genotype))
        summary=summary[summary > 2]
        if(length(unique(draw$Genotype)) == 1 | length(summary)==1 ){
            return(NULL)
        }else{            
            if(opt$btl==1){
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
            compare=combn(names(summary),m=2,simplify=F)
            p=ggboxplot(draw,x="Genotype",y="Pheno",group="Genotype",fill="Genotype")+stat_compare_means(comparisons=compare,aes(label=..p.format..),method="t.test")+scale_fill_npg()
            p=p+theme_bw()
            ggsave(file=paste(phe.name,result$peak[x],"pxg.pdf",sep="."),p)
            ggsave(file=paste(phe.name,result$peak[x],"pxg.png",sep="."),p)
            df=compare_means(data=draw[draw$Genotype %in% unlist(compare),],Pheno~Genotype,method="t.test")
            colnames(df)[1]<-"Pheno"
            df[,1] <- phe.name
            avg <- draw %>% group_by(Genotype) %>% summarise(mean=mean(Pheno),sd=sd(Pheno))
            avgv <- avg$mean
            sd <- avg$sd
            names(avgv) <- avg$Genotype
            names(sd) <- avg$Genotype
            df$group1num=summary[df$group1]
            df$group1avg=avgv[df$group1]
            df$group1sd=sd[df$group1]
            df$group2num=summary[df$group2]
            df$group2avg=avgv[df$group2]
            df$group2sd=sd[df$group2]
            df$peak=result$peak[x]
            df
            }
        }
    }
    compare=lapply(1:nrow(result),drawpxg,result=result,d=d,phe.name=phe.name);
    compare_df=do.call(rbind,compare)
    write.table(file=paste(phe.name,"peaks.pxg.result",sep="."),compare_df,sep="\t",row.names=F)
}else{
    cat("",file=paste(phe.name,"qtl-result.result",sep="."))
}

    print("haha2")
    # 将chr字段提取为数字进行排序
    scan$chr_num <- as.numeric(gsub("chr", "", scan$chr))
    # 排序时使用chr_num进行排序
    scan = scan[order(scan$chr_num, scan$pos), ]
    # 清除临时生成的chr_num字段
    scan$chr_num <- NULL
    legend = round(pm.result, 2)
    # 绘图部分
    pdf(file=paste(phe.name,"scan.pdf",sep="."), height=9, width=16)
    plot(scan)
    print("haha")
    head(pm.result)
    abline(h=pm.result, col=rainbow(length(pm.result)))
    legend("topright", legend=legend, col=rainbow(length(pm.result)))
    dev.off()
    png(file=paste(phe.name,"scan.png",sep="."), height=900, width=1600)
    plot(scan)
    abline(h=pm.result, col=rainbow(length(pm.result)))
    legend("topright", legend=legend, col=rainbow(length(pm.result)))
    dev.off()

    # scan=scan[order(as.numeric(as.character(scan$chr)),as.numeric(scan$pos)),]
    # legend=round(pm.result,2)
    # pdf(file=paste(phe.name,".scan.pdf",sep=""),height=9,width=16)
    # plot(scan)
    # abline(h=pm.result,col=rainbow(length(pm.result)))
    # legend("topright",legend=legend,col=rainbow(length(pm.result)))
    # dev.off()
    # png(file=paste(phe.name,".scan.png",sep=""),height=900,width=1600)
    # plot(scan)
    # abline(h=pm.result,col=rainbow(length(pm.result)))
    # legend("topright",legend=legend,col=rainbow(length(pm.result)))
    # dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
