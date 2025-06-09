#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'bulksize','b',0,'character',
	'bootstrap','t',0,'character',
	'popstruc','p',0,'character',
	'qvalue','q',0,'character',
    'mutmap','m',0,'character',
    'noboot','n',0,'logical',
	'help','h',0,'logical',
    'abs','a',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--infile	the sliding win  file
	--outfile	the output file
	--bulksize	the bulksize  default 30
	--bootstrap	the bootstrap numbers
	--popstruc	the population structure default F2
	--qvalue	the qvalue of quantile default 0.05
	--abs	the index is abs or not
    --noboot	use quantile instead of bootstrap
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }
if ( is.null(opt$bulksize)){opt$bulksize=30}
if ( is.null(opt$seqdep)){opt$seqdep=30}
if ( is.null(opt$popstruc)){opt$popstruc="F2"}
if ( is.null(opt$qvalue)){ opt$qvalue=0.005}
if ( is.null(opt$bootstrap)){ opt$bootstrap=1000}
if ( is.null(opt$noboot)){ opt$noboot=F}
if(is.null(opt$abs)){opt$abs=F}
if(is.null(opt$mutmap)){opt$mutmap=F}
opt$bootstrap=as.numeric(opt$bootstrap)
opt$qvalue=1-as.numeric(opt$qvalue)
library(dplyr)
simulateAlleleFreq <- function(n, pop = "F2") {
    if (pop == "F2") {
        mean(sample(
            x = c(0, 0.5, 1),
            size = n,
            prob = c(1, 2, 1),
            replace = TRUE
        ))
    } else if (pop == "F1"){
            mean(sample(x = c(0, 0.5, 1),
            size = n,
            prob = c(1,2,1),
            replace = TRUE
    ))
    }else {
        mean(sample(
            x = c(0, 1),
            size = n,
            prob = c(1, 1),
            replace = TRUE
        ))
    }
}
simulateSNPindex <-
    function(depth,
             altFreq1,
             altFreq2,
             replicates = 1000,
             filter = NULL,abs=F,mutmap=F) {

        SNPindex_H <- rbinom(replicates, size = depth, altFreq1) / depth
        SNPindex_L <- rbinom(replicates, size = depth, altFreq2) / depth
        deltaSNP <- SNPindex_H - SNPindex_L

        if (!is.null(filter)) {
            deltaSNP <- deltaSNP[SNPindex_H >= filter | SNPindex_L >= filter]
        }
	if (mutmap){
		SNPindex_H
	}else{
		if(abs){
			abs(deltaSNP)
		}else{
			deltaSNP
		}
	}
    }
simulateConfInt <- function(popStruc = "F2",
                            bulkSize,
                            depth = 1:100,
                            replications = 10000,
                            filter = NULL,
                            intervals = seq(0, 1,0.0001),
			    abs=NULL,
			    mutmap=NULL
			    ) {
    if (popStruc == "F2") {
        message(
            "Assuming bulks selected from F2 population, with ",
            paste(bulkSize, collapse = " and "),
            " individuals per bulk."
        )
    } else {
        message(
            "Assuming bulks selected from F1 population, with ",
            bulkSize,
            " individuals per bulk."
        )
    }

    if (length(bulkSize) == 1) {
        message("The 'bulkSize' argument is of length 1, setting number of individuals in both bulks to: ", bulkSize)
        bulkSize[2] <- bulkSize[1]
    }

    if (length(bulkSize) > 2) {
        message("The 'bulkSize' argument is larger than 2. Using the first two values as the bulk size.")
    }

    if (any(bulkSize < 0)) {
        stop("Negative bulkSize values")
    }

    #makes a vector of possible alt allele frequencies once. this is then sampled for each replicate
    print(replications);
    tmp_freq <-
        replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize[1], pop = popStruc))

    tmp_freq2 <-
        replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize[2], pop = popStruc))

    message(
        paste0(
            "Simulating ",
            replications,
            " SNPs with reads at each depth: ",
            min(depth),
            "-",
            max(depth)
        )
    )
    CI <- sapply(
        X = depth,
        FUN = function(x)
        {
            quantile(
                na.rm = TRUE,
                x = simulateSNPindex(
                    depth = ceiling(x),
                    altFreq1 = sample(
                        x = tmp_freq,
                        size = replications,
                        replace = TRUE
                    ),
                    altFreq2 = sample(
                        x = tmp_freq2,
                        size = replications,
                        replace = TRUE
                    ),
                    replicates = replications,
                    filter =filter,
		    abs=abs,
		    mutmap=mutmap
                ),
                probs = intervals,
                names = TRUE
            )
        }
    )
    CI <- as.data.frame(CI)

    CI <- data.frame(t(CI))
    print("haha")

    print(CI)
    names(CI) <- intervals
    CI <- cbind(depth, CI)

    #to long format for easy plottingCI$
    # tidyr::gather(data = CI,
    #     key = interval,
    #     convert = TRUE,
    #     value = SNPindex,-depth) %>%
    #     dplyr::mutate(Confidence = factor(ifelse(
    #         interval > 0.5,
    #         paste0(round((1 - interval) * 200, digits = 1), "%"),
    #         paste0((interval * 200), "%")
    # )))
    CI
}


opt$bulksize=as.numeric(opt$bulksize)
if(opt$popstruc=="F1"){opt$abs=T}
data<-read.table(opt$infile,head=TRUE,comment.char="^")
depth=as.numeric(levels(as.factor(data$mdepth)))
opt$qvalue= as.numeric(opt$qvalue)
CI<- simulateConfInt(
        popStruc = opt$popstruc,
        bulkSize = opt$bulksize,
        depth = depth,
        replications = opt$bootstrap,
        intervals = seq(0+1/opt$bootstrap,1,1/opt$bootstrap),
	    abs=opt$abs,
	    mutmap=opt$mutmap
)
write.table(file=paste(opt$outfile,"bootstrap.detail",sep="."),CI,quote=F,row.names=F,sep="\t")
print(colnames(CI))
if(opt$noboot){
    qt <- quantile(data$slidingD, opt$qvalue, na.rm = T)
    nCI<-data.frame(depth=CI$depth, CI=qt)
}else{
    nCI<-data.frame(depth=CI$depth, CI=CI[colnames(CI)==opt$qvalue])
}
colnames(nCI)=c("mdepth","CI")
ndata <- left_join(data,nCI,by="mdepth")
print(colnames(ndata))
ndata <- arrange(ndata,X.chr,pos1,pos2)
print(opt$mutmap)
if (!("ncount" %in% colnames(ndata))){
    ndata$ncount=1;
}
if(opt$mutmap){
outd=ndata %>% select(X.chr,pos,delta,pos1,pos2,slidingD,ncount,CI)
}else{
outd=ndata %>% select(X.chr,pos,index1,index2,delta,pos1,pos2,slidingI1,slidingI2,slidingD,ncount,CI)
}
write.table(file=paste(opt$outfile,"bootstrap.result",sep="."),outd,quote=F,row.names=F,sep="\t")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
