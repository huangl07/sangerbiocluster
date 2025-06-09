#!/usr/bin/env Rscript
times <- Sys.time()
library("getopt")
library(plotrix)
options(bitmapType = "cairo")
spec <- matrix(c(
    "i", "i", 1, "character",
    "o", "o", 1, "character",
    "s", "s", 2, "character",
    "help", "c", 0, "logical",
    "d", "d", 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec = NULL) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("
Usage example:
	Rscript base_distrbution.r --i  --o

Usage:
	--i     base_distrbution
	--o	     base_distrbution picture name
	--s	single
    --d sample name
	--help		usage
\n")
    q(status = 1)
}

if (!is.null(opt$help)) {
    print_usage(spec)
}
if (is.null(opt$i)) {
    print_usage(spec)
}
if (is.null(opt$o)) {
    print_usage(spec)
}

xmax <- 2000
a <- read.table(opt$i)
total_counts <- sum(a[, 2])
cum_counts <- cumsum(a[, 2])
sum_depth <- sum(a[, 1] * a[, 2])
median_pos <- (total_counts + 1) / 2
ave <- sum_depth / total_counts
if (total_counts %% 2 != 1) {
    median_value <- a[which.min(abs(cum_counts - median_pos)), 1]
} else {
    # 如果数据集的数量是偶数，则取中间两个数值的平均值
    lower_value <- a[which.min(abs(cum_counts - median_pos)), 1]
    upper_value <- a[which.min(abs(cum_counts - (median_pos - 1))), 1]
    median_value <- (lower_value + upper_value) / 2
}
a <- a[-nrow(a), ]
b <- cumsum(as.numeric(a[, 2]))
sums <- sum(as.numeric(a[, 2]))
b <- b / sums * 100
c <- a[, 2] / sums * 100
write.table(
    data.frame(
        sample = opt$d,
        mean = ave,
        median = median_value
    ),
    paste(opt$o, ".median", sep = ""),
    quote = F, row.names = F, col.names = T, sep = "\t"
)

if (is.null(opt$s)) {
    lmax <- max(c)
    xmean <- a$V1[match(lmax, c)]
    if (xmean < 10) {
        xmax <- 20
    } else if (xmean < 100) {
        xmax <- 200
    } else {
        xmax <- 1000
    }
    pdf(paste(opt$o, ".pdf", sep = ""))
    twoord.plot(a[, 1], c, a[, 1], b, xlim = c(0, xmax), lylim = c(0, lmax * 1.2), rylim = c(0, 110), lcol = "red3", rcol = "blue", xlab = "Sequencing depth", ylab = "Percent of base", rylab = "Percent of cumulative base", type = c("l", "l"))
    legend("right", legend = c("Percent of base (%)", "Percent of cumulative base (%)"), col = c("red3", "blue"), lty = 1, bty = "n", text.col = c("red3", "blue"), cex = 0.5)
    title(main = paste("Depth of", opt$d))
    dev.off()
    png(paste(opt$o, ".png", sep = ""))
    twoord.plot(a[, 1], c, a[, 1], b, xlim = c(0, xmax), lylim = c(0, lmax * 1.2), rylim = c(0, 110), lcol = "red3", rcol = "blue", xlab = "Sequencing depth", ylab = "Percent of base", rylab = "Percent of cumulative base", type = c("l", "l"))
    legend("right", legend = c("Percent of base (%)", "Percent of cumulative base (%)"), col = c("red3", "blue"), lty = 1, bty = "n", text.col = c("red3", "blue"), cex = 0.5)
    title(main = paste("Depth of", opt$d))
    dev.off()
} else {
    lmax <- max(c[-c(1:5)])
    xmean <- a$V1[match(lmax, c)]
    if (xmean < 10) {
        xmax <- 50
    } else if (xmean < 100) {
        xmax <- 500
    } else {
        xmax <- 1000
    }
    pdf(paste(opt$o, ".pdf", sep = ""))
    plot(a[, 1], c, col = "blue", type = "l", xlab = "Sequencing depth", ylab = "Percent of base (%)", xlim = c(0, xmax))
    title(main = paste("Depth of", opt$d))
    dev.off()
    png(paste(opt$o, ".png", sep = ""))
    plot(a[, 1], c, col = "blue", type = "l", xlab = "Sequencing depth", ylab = "Percent of base (%)", xlim = c(0, xmax))
    title(main = paste("Depth of", opt$d))
    dev.off()
}
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
