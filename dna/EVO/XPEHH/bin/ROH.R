library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-p", "--pop"), type = "character", action = "store", default = NULL,help = "Input keys word by pop"),
  make_option(c("-o", "--threshold"), type = "integer", action = "store",  default = 0.7 , help = "output file")

)
opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"))
if(is.null(opt$pop)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$threshold)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
times<-Sys.time()

library(detectRUNS)
library(dplyr)
library(ggplot2)
slidingRuns<-slidingRUNS.run(genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."))
summaryList <- summaryRuns(runs = slidingRuns,genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."))
table<-tableRuns(slidingRuns,genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."),threshold=opt$threshold)
mapFile=paste(opt$pop,"map",sep=".")
genotypeFile=paste(opt$pop,"ped",sep=".")
runs=slidingRuns
names(runs) <- c("POPULATION", "IND", "CHROMOSOME", "COUNT","START", "END", "LENGTH")
BP <- NULL
P <- NULL
CHR <- NULL
if (file.exists(mapFile)) {
    mappa <- data.table::fread(mapFile, header = F)
}else {
    stop(paste("file", mapFile, "doesn't exists"))
}
names(mappa) <- c("CHR", "SNP_NAME", "x", "POSITION")
mappa$x <- NULL
print("Calculation % SNP in ROH")
all_SNPinROH <- data.frame(SNP_NAME = character(), CHR = integer(), 
    POSITION = numeric(), COUNT = integer(), BREED = factor(), 
    PERCENTAGE = numeric(), stringsAsFactors = FALSE)
total <- length(unique(runs$CHROMOSOME))
print(paste("Chromosome founds: ", total))
n = 0
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (chrom in sort(unique(runs$CHROMOSOME))) {
    runsChrom <- runs[runs$CHROMOSOME == chrom, ]
    mapChrom <- mappa[mappa$CHR == chrom, ]
    snpInRuns <- detectRUNS:::snpInsideRunsCpp(runsChrom, mapChrom, genotypeFile)
    all_SNPinROH <- rbind.data.frame(all_SNPinROH, snpInRuns)
    n = n + 1
    setTxtProgressBar(pb, n)
}
close(pb)
print("Calculation % SNP in ROH finish")
print("Manhattan plot: START")
group_list = unique(all_SNPinROH$BREED)

    for (group in group_list) {
        print(paste("Processing Groups:", group))
        subset_group = subset(all_SNPinROH, all_SNPinROH$BREED == 
            group)
        names(subset_group) <- c("SNP", "CHR", "BP", "COUNT", 
            "GROUP", "P")
        dfpos = subset_group
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

        axisdf <- dpos %>% group_by(CHR) %>% summarize(center=(as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )
        chrNum=nlevels(as.factor(dpos$CHR))
        p <- ggplot(dpos)
        p <- p + geom_point(aes(x = BPcum, y = P, colour = as.factor(CHR)), 
            alpha = 2/3)
        p <- p + scale_color_manual(values = rep(c("red", "blue"), 
            round(chrNum/2, 0) + 1))
        p <- p + scale_size(range = c(0.1, 0.1)) + ylim(0, 100)
        p <- p + theme_bw(base_size = 11) + theme(legend.position = "none")
        p <- p + scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center )
        roh_plot <- p + ggtitle("manhattan Plot") + xlab("Chromosome") + 
            ylab("% SNP in Runs") + theme(plot.title = element_text(hjust = 0.5))
        ggsave(filename = paste("ROH_", group, 
                ".pdf", sep = ""), plot = roh_plot, device = "pdf",width=16,height=9)
        ggsave(filename = paste("ROH_", group, 
                ".png", sep = ""), plot = roh_plot, device = "png",width=16,height=9)
    }
write.table(file="ROH.sliding",slidingRuns)
write.table(file="ROH.table",table)
write.table(file="summary_ROH_count_chr",summaryList$summary_ROH_count_chr)
write.table(file="summary_ROH_percentage_chr",summaryList$summary_ROH_percentage_chr)
write.table(file="summary_ROH_count",summaryList$summary_ROH_count)
write.table(file="summary_ROH_percentage",summaryList$summary_ROH_percentage)
write.table(file="summary_ROH_mean_chr",summaryList$summary_ROH_mean_chr)
write.table(file="summary_ROH_mean_class",summaryList$summary_ROH_mean_class)
write.table(file="result_Froh_genome_wide",summaryList$result_Froh_genome_wide)
write.table(file="result_Froh_chromosome_wide",summaryList$result_Froh_chromosome_wide)
write.table(file="result_Froh_class",summaryList$result_Froh_class)
write.table(file="SNPinRun",summaryList$SNPinRun)


escaptime=Sys.time()-times
print("Done!")
print(escaptime)
