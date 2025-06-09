library(lme4)
library(lsmeans)
blup<-function(value){
    v=csv[[value]]
    g=csv[["trt"]]
    r=csv[["replicates"]]
    m1=lmer(v~g+(1|r))
    re=emmeans(m1,"g")
    re=as.data.frame(re)
    df=data.frame(id=re$g)
    df[[value]]=re$emmean
    df
}

csv<-read.csv("YHJ.csv")
list=lapply(colnames(csv)[3:ncol(csv)],blup)
do.call(list,left_join(by="id"))

