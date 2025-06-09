library(vegan)
library(lfmm)
library(qvalue)
library(adegenet)

gen<-read.PLINK("pop.raw",mapfile="pop.map")
gen.imp <- apply(gen, 2,
    function(x){
        replace(x,
        is.na(x),
        as.numeric(names(which.max(table(x))))
    )
    }
)
env<-read.csv("env.txt",sep="\t")
env[,1]=as.character(env[,1])

senv=env[env$individual %in% rownames(gen.imp),]
gen.imp.re=gen.imp[senv$individual,]
senv[,1]=as.character(senv[,1])
gen.imp=gen.imp.re[senv[,1],]
env=senv
pred<-env[,3:9]
pred.pca <- rda(pred, scale=T)
png("ScreenplotENV.png")
screeplot(pred.pca)
dev.off()
gen.pca <- rda(gen.imp, scale=T)
png("ScreenplotGene.png")
screeplot(gen.pca)
dev.off()


pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)
myfun<-function(n){
    K=n;
    lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=K)
    pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=lfmm, calibrate="gif")
    png(paste("lfmm","K",K,"png",sep="."))
    hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")
    dev.off()
    pv$gif
}
result=data.frame(K=1:6)
result$gif=sapply(1:6,myfun)
bestK=result$K[which.min(abs(result$gif-1))]
lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=bestK)
pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=lfmm, calibrate="gif")
qv <- qvalue(pv$calibrated.pvalue)$qvalues
qv=as.data.frame(qv)
myfun1<-function(n){
    pred.PC1<-senv[,n]
    result=data.frame(K=1:6)
    result$gif=sapply(1:6,myfun)
    bestK=result$K[which.min(abs(result$gif-1))]  
    lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=bestK)
    pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=lfmm, calibrate="gif")
    qvalue(pv$calibrated.pvalue)$qvalues[,1]
}
qv[[colnames(senv)[3]]]=myfun1(3)
qv[[colnames(senv)[4]]]=myfun1(4)
qv[[colnames(senv)[5]]]=myfun1(5)
qv[[colnames(senv)[6]]]=myfun1(6)
qv[[colnames(senv)[7]]]=myfun1(7)
qv[[colnames(senv)[8]]]=myfun1(8)
qv[[colnames(senv)[9]]]=myfun1(9)

