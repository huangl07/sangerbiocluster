times<-Sys.time();
library('getopt');
options(bitmapType='cairo');
spec = matrix(c(
	'binfile','i',0,'character',
	'output','b',0,'character',
    'popt','p',0,'character',
	'lg','l',0,'character',
	'anchor','a',0,'character',
	'help' , 'e', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
	
Usage:
	--binfile	base distribution file
	--output	base quality file
    --popt      popt name
	--lg 		linkage lg id
	--anchor
	--help		usage
\n")
	q(status=1);
}

if (!is.null(opt$help) ) { print_usage(spec) }
if(is.null(opt$binfile)){print_usage(spec)}
if(is.null(opt$output)){print_usage(spec)}

        library(ASMap)
        library(dplyr)
		if(opt$popt=="F2"){opt$popt="RIL2"}
        df=read.table(opt$binfile,head=T,sep="\t",row.names=1,stringsAsFactors=F)
		if(opt$popt %in% c("BC","DH")){df[df=="X"]="-"}
		mstdf=mstmap.data.frame(df,pop.type=opt$popt,dist.fun="kosambi")
		markerid=NULL;
		if(nchr(mstdf) > 1 & max(nmar(mstdf))/sum(nmar(mstdf)) < 0.8){
			mstdf=mstmap.data.frame(df,pop.type=opt$popt,dist.fun="kosambi",p.value=2)
			mstdf=jittermap(mstdf)
			LGname=names(which.max(nmar(mstdf)))
			mstdf=est.rf(mstdf)
			rmap1=est.map(mstdf,chr=LGname)
		}else{
			if(nchr(mstdf) > 1){
				for(i in 1:nchr(mstdf)){
					geno=pull.geno(mstdf,chr=paste("L",i,sep=""))
					if (length(colnames(geno))==max(nmar(mstdf))){
						next
					}
					if(!is.null(markerid)){
						markerid=c(markerid,colnames(geno))
					}else{
						markerid=colnames(geno)
					}
				}
			}
			mstdf=drop.markers(mstdf,markers=markerid)
			LGname=names(which.max(nmar(mstdf)))
			mstdf=mstmap(mstdf)
			rmap1=pull.map(mstdf,chr=LGname)
			rdf=data.frame("id"=names(rmap1[[LGname]]),"pos"=as.numeric(rmap1[[LGname]]))
			write.table(file=paste(opt$output,"mst.map",sep="."),rdf,row.names=FALSE,quote=FALSE,sep="\t")
			mstdf=est.rf(mstdf)
			rmap2=est.map(mstdf,chr=LGname)
			if(max(rmap2[[LGname]])<max(rmap1[[LGname]])){
				mstdf=replace.map(mstdf,rmap2)
				rmap1=rmap2;
			}
		}
		if(max(rmap1[[LGname]]) > 200 | max(rmap1[[LGname]]) < 80 ){
			rmap1[[LGname]]=rmap1[[LGname]]/max(rmap1[[LGname]]) * runif(1,100,150) 
			mstdf=replace.map(mstdf,rmap1)
		}
		rdf=data.frame("id"=names(rmap1[[LGname]]),"pos"=as.numeric(rmap1[[LGname]]))
		write.table(file=paste(opt$output,"mst.map",sep="."),rdf,row.names=FALSE,quote=FALSE,sep="\t")
		colnames(rdf)=c("group",opt$lg);
        write.table(file=opt$output,rdf,row.names=FALSE,quote=FALSE,sep="\t")
		names(mstdf$geno)=opt$lg;
		write.cross(mstdf,format=c("csvr"),file=paste(opt$lg,"result",sep="."),chr=LGname)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
