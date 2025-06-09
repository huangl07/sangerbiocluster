library(purrr)
library(dndscv)
a<-read.table("CT.list")
ndf=map_dfr(a$V2,read.table,head=F,sep="\t",quote="",.id="id")
ndfs=ndf[,c(1,2,3,4,5)]
ndfs=ndfs[!duplicated(ndfs),]
ndfs$id=a$V1[as.numeric(ndfs$id)]

dndsout = dndscv(ndfs,refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Hdndsout = dndscv(ndfs[grepl("H",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Ldndsout = dndscv(ndfs[grepl("L",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Zdndsout = dndscv(ndfs[grepl("Z",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Mdndsout = dndscv(ndfs[grepl("M",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)

write.table(file="CT_pan.dndsout.result",dndsout$sel_cv)
write.table(file="CT_H.dndsout.result",Hdndsout$sel_cv)
write.table(file="CT_L.dndsout.result",Ldndsout$sel_cv)
write.table(file="CT_Z.dndsout.result",Zdndsout$sel_cv)
write.table(file="CT_M.dndsout.result",Mdndsout$sel_cv)


a<-read.table("C.list")
ndf=map_dfr(a$V2,read.table,head=F,sep="\t",quote="",.id="id")
ndfs=ndf[,c(1,2,3,4,5)]
ndfs=ndfs[!duplicated(ndfs),]
ndfs$id=a$V1[as.numeric(ndfs$id)]

dndsout = dndscv(ndfs,refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Hdndsout = dndscv(ndfs[grepl("H",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Ldndsout = dndscv(ndfs[grepl("L",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Zdndsout = dndscv(ndfs[grepl("Z",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Mdndsout = dndscv(ndfs[grepl("M",ndfs$id),],refdb="~/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)

write.table(file="C_pan.dndsout.result",dndsout$sel_cv)
write.table(file="C_H.dndsout.result",Hdndsout$sel_cv)
write.table(file="C_L.dndsout.result",Ldndsout$sel_cv)
write.table(file="C_Z.dndsout.result",Zdndsout$sel_cv)
write.table(file="C_M.dndsout.result",Mdndsout$sel_cv)


 myfun<-function(dndsout,out){
    sel_cv=dndsout$sel_cv
    dir.create(out)
    draw= sel_cv[sel_cv$qmis_cv <0.05 |sel_cv$qtrunc_cv <0.05 ,]
    draw$type[draw$qmis_cv < 0.05 & draw$qtrunc_cv > 0.05]="missense"
    draw$type[draw$qmis_cv > 0.05 & draw$qtrunc_cv < 0.05]="truncating"
    draw$type[draw$qmis_cv < 0.05 & draw$qtrunc_cv < 0.05]="both"
    p1=ggplot(draw)+geom_point(aes(x=wnon_cv,y=wmis_cv,col=type))+xlab("dN/dS(truncating)")+ylab("dN/dS(missense)")+theme_bw()
    ggsave(p1,file=paste(out,"fig2-B.pdf",sep="/"))
    ggsave(p1,file=paste(out,"fig2-B.png",sep="/"))
    all=sel_cv[ sel_cv$gene_name %in% dndsout$annotmuts$gene,]
    all$wmis_cv[all$wmis_cv > 10]=10
    p2=ggplot(all)+geom_histogram(aes(x=wmis_cv),binwidth=0.1)+ylab("Genes")+xlab("dN/dS(missense)")+xlim(c(0,10.05))
    ggsave(p2,file=paste(out,"fig3-A.pdf",sep="/"))
    ggsave(p2,file=paste(out,"fig3-A.png",sep="/"))
    model =glm.nb(n_syn~offset(log(exp_syn)),data=dndsout$genemuts)
    pred=data.frame(geneid=dndsout$genemuts$gene_name,predict=predict(model,newdata=dndsout$genemuts))
    p=ggplot(pred)+geom_histogram(aes(x=predict),binwidth=0.1)+xlim(c(0,5))
    ggsave(file=paste(out,"Figure3-B.pdf",sep="/"),p)
    ggsave(file=paste(out,"Figure3-B.png",sep="/"),p)
    draw$group="NONE"
    draw$group[draw$wmis_cv <= 0.25]="<0.25"
    draw$group[draw$wmis_cv <= 0.75 & draw$wmis_cv > 0.25]="<0.75"
    draw$group[draw$wmis_cv >= 1.5 & draw$wmis_cv < 5]=">1.5"
    draw$group[draw$wmis_cv >=5]=">5"
    p3=ggplot(draw[draw$group != "NONE",])+geom_bar(aes(x=group,fill=group),stat="count")+theme_bw()+ylab("% of genes")+xlab("dN/DS")
    ggsave(p3,file=paste(out,"fig3-C.pdf",sep="/"))
    ggsave(p3,file=paste(out,"fig3-C.png",sep="/"))
    p3=ggplot(draw[draw$group != "NONE" & draw$gene_name %in% all$gene_name,])+geom_bar(aes(x=group,fill=group),stat="count")+theme_bw()+ylab("% of genes")+xlab("dN/DS")
    ggsave(p3,file=paste(out,"fig3-C-have-variant.pdf",sep="/"))
    ggsave(p3,file=paste(out,"fig3-C-have-variant.png",sep="/"))
    p3=ggplot(draw[draw$group != "NONE" & draw$gene_name %in% all$gene_name & draw$wmis_cv != 0,])+geom_bar(aes(x=group,fill=group),stat="count")+theme_bw()+ylab("% of genes")+xlab("dN/DS")
    ggsave(p3,file=paste(out,"fig3-C-have-variant-no0.pdf",sep="/"))
    ggsave(p3,file=paste(out,"fig3-C-have-variant-no0.png",sep="/"))
}

myfun(dndsout,"A")
myfun(Hdndsout,"H")
myfun(Zdndsout,"Z")
myfun(Ldndsout,"L")
myfun(Mdndsout,"M")

sel_cv=dndsout$sel_cv
Hsel_cv=Hdndsout$sel_cv
Zsel_cv=Zdndsout$sel_cv
Lsel_cv=Ldndsout$sel_cv
Msel_cv=Mdndsout$sel_cv

sel_cv$group="pan"
Hsel_cv$group="H"
Zsel_cv$group="Z"
Lsel_cv$group="L"
Msel_cv$group="M"

genemuts=dndsout$genemuts
Hgenemuts=Hdndsout$genemuts
Zgenemuts=Zdndsout$genemuts
Lgenemuts=Ldndsout$genemuts
Mgenemuts=Mdndsout$genemuts
genemuts$group="pan"
Hgenemuts$group="H"
Zgenemuts$group="Z"
Lgenemuts$group="L"
Mgenemuts$group="M"

draw=rbind(genemuts,Hgenemuts,Zgenemuts,Lgenemuts,Mgenemuts)

draw$global=(draw$n_non+draw$n_spl+draw$n_miss)/(draw$exp_miss+draw$exp_spl+draw$exp_non)

draw=draw[draw$gene_name %in% dndsout$annotmuts$gene,]
p=ggplot(draw[draw$global != 0,])+geom_boxplot(aes(x=group,y=global))
ggsave("figure-4-B-no0.pdf",p)

gene<-read.table("~/sg-users/long.huang/DNDSTOOLS/mouse.gene",head=F,sep="\t",)
geneid=paste(gene$V1,gene$V2,sep=":")
draw=rbind(genemuts,Hgenemuts,Zgenemuts,Lgenemuts,Mgenemuts)
draw$global=(draw$n_non+draw$n_spl+draw$n_miss)/(draw$exp_miss+draw$exp_spl+draw$exp_non)
draw=draw[draw$gene_name %in% geneid,]
p=ggplot(draw[draw$global != 0,])+geom_boxplot(aes(x=group,y=global))
ggsave("figure-4-A-no0.pdf",p)


library(purrr)
library(dndscv)
library(ggplot2)
library(MASS)
a<-read.table("../CT26.xls")
ndf=map_dfr(a$V2,read.table,head=T,sep="\t",quote="",.id="id")
ndfs=ndf[,c(1,2,3,4,5)]
ndfs=ndfs[!duplicated(ndfs),]
ndfs$id=a$V1[as.numeric(ndfs$id)]
dndsout = dndscv(ndfs,refdb="/mnt/lustre/users/sanger-dev/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Hdndsout = dndscv(ndfs[grepl("H",ndfs$id),],refdb="/mnt/lustre/users/sanger-dev/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Ldndsout = dndscv(ndfs[grepl("L",ndfs$id),],refdb="/mnt/lustre/users/sanger-dev/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Zdndsout = dndscv(ndfs[grepl("Z",ndfs$id),],refdb="/mnt/lustre/users/sanger-dev/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
Mdndsout = dndscv(ndfs[grepl("M",ndfs$id),],refdb="/mnt/lustre/users/sanger-dev/sg-users/long.huang/DNDSTOOLS/example_output_refcds.rda", cv=NULL,outmat=T,constrain_wnon_wspl=T,outp=3,max_coding_muts_per_sample=5000)
sel_cv=dndsout$sel_cv
Hsel_cv=Hdndsout$sel_cv
Zsel_cv=Zdndsout$sel_cv
Lsel_cv=Ldndsout$sel_cv
Msel_cv=Mdndsout$sel_cv
sel_cv$group="pan"
Hsel_cv$group="H"
Zsel_cv$group="Z"
Lsel_cv$group="L"
Msel_cv$group="M"
write.table(sel_cv, "pan.sel_cv.xls",sep="\t",row.names=F,quote=F)
write.table(Hsel_cv, "H.sel_cv.xls",sep="\t",row.names=F,quote=F)
write.table(Zsel_cv, "Z.sel_cv.xls",sep="\t",row.names=F,quote=F)
write.table(Lsel_cv, "L.sel_cv.xls",sep="\t",row.names=F,quote=F)
write.table(Msel_cv, "M.sel_cv.xls",sep="\t",row.names=F,quote=F)