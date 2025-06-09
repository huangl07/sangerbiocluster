library(optparse)
options(bitmapType = "cairo")
option_list <- list(
    make_option(c("-x", "--xpehh"),
        type = "character",
        action = "store", help = "Input xpclr file"
    ),
    make_option(c("-t", "--threshold"),
        type = "numeric",
        action = "store",
        default = 0.05, help = "threshold for select"
    )
)
opt <- parse_args(OptionParser(
    option_list = option_list,
    usage = "Usage: %prog [options] \nDescription: This Script is used to draw xpclr manhattan!"
))
if (is.null(opt$xpehh)) {
    opt <- parse_args(OptionParser(
        option_list = option_list,
        usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"
    ), args = "--help")
}

times <- Sys.time()

library(ggplot2)
library(dplyr)
library(ggpubr)


dfs <- read.table(opt$xpehh, sep = "\t", header = TRUE)
#         sed -i  '1i Chromosome\\tStart\\tEnd\\tNumSites\\tFrac\\tPecentile\\tiHS" total.ihs.windows.result
### sed '1i Compare\\tChromosome\\tPosition1\\tPosition2\\tLikelihood\\tAlpha\\tStartPos\\tEndPos'> total.clr.result


colnames(dfs) <- c("compare","id", "pos1", "pos2", "likelihood","alpha","StartPos","EndPos")
dfs <- dfs[grep("chr", dfs$id), ]
df1 <- dfs %>% filter(likelihood > quantile(as.numeric(dfs$likelihood), 1 - opt$threshold, na.rm = TRUE))
pop <- gsub(".likelihood.windows.result", "", opt$xpehh)
print(pop)
write.table(file = paste(pop, "likelihood.select", sep = "."), df1, row.names = FALSE, quote = FALSE, sep = "\t")
draw=function(compare){
    df <- dfs[compare == compare]
    print(compare)
    df$pos <- (df$pos1 + df$pos2) / 2
    df$BP <- df$pos
    print(nrow(df))
    df$likelihood <- as.numeric(df$likelihood)
    quantile1 <- quantile(df$likelihood, 1 - opt$threshold, na.rm = TRUE)
    quantile2 <- quantile(df$likelihood, opt$threshold, na.rm = TRUE)
    print("haha")
    print(quantile1)
    print(quantile2)
    df$CHR <- df$id
    lev <- NULL
    lev$CHR <- levels(as.factor(df$CHR))
    lev$order <- stringr::str_rank(lev$CHR, numeric = TRUE)
    dfpos <- merge(df, lev, by = "CHR")
    dfpos <- arrange(dfpos, order, BP)
    print("haha")
    dpos <- dfpos %>%
        group_by(order) %>%
        summarise(chr_len = 1.2 * max(BP)) %>%
        mutate(tot = cumsum(chr_len) - chr_len) %>%
        select(-chr_len) %>%
        left_join(dfpos, ., by = c("order" = "order")) %>%
        arrange(order, BP) %>%
        mutate(BPcum = BP + tot)
    axisdf <- dpos %>%
        group_by(CHR) %>%
        summarize(center = (as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2)
    dpos$col <- 1
    dpos$col[dpos$order %% 2 == 0] <- 2
    dpos$col[dpos$likelihood > quantile1] <- 3
    p1 <- ggplot(dpos) +
        geom_point(aes(x = BPcum, y = likelihood, color = as.factor(col))) +
        geom_hline(yintercept = quantile1, linetype = "dashed", col = "black")
    p1 <- p1 + scale_color_manual(values = rep(c("slategray3", "skyblue4", "red"), nlevels(as.factor(dpos$CHR)))) +
        scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
        ylab("CLR") + labs(title = paste("CLR along All Chromosome")) + xlab("Chromosome") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.y.left = element_line(),
            axis.line.x.bottom = element_line(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        )

    ggsave(p1, file = paste(compare, "clr.manhattan.png", sep = "."), device = "png", height = 9, width = 16)
    ggsave(p1, file = paste(compare, "clr.manhattan.pdf", sep = "."), device = "pdf", height = 9, width = 16)
}
levels=levels(factor(dfs$compare))
print(levels)
lapply(levels,draw)

escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
