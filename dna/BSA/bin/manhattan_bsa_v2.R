#!/usr/bin/env Rscript

times <- Sys.time()

# load library
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggh4x"))
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("latex2exp"))
options(bitmapType = "cairo")
#########################################################################################
options(scipen = 200)
spec <- matrix(c(
	'result', 'r', 1, 'character',
	'threshold', 't', 0, 'numeric',
	'chr', 'c', 1, 'character',
	'output', 'o', 0, 'character',
	'pcol', 'p', 0, 'character',
	'lcol', 'l', 0, 'character',
	'qcol', 'q', 0, 'character',
	'xlab', 'x', 0, 'character',
	'ylab', 'y', 1, 'character',
	'chrgap', 'g', 0, 'numeric',
      'mutmap', 'm', 0, 'logical',
    'abs','a',0,'logical',
	'help', 'h', 0, 'logical'), byrow = T, ncol = 4)
opt <- getopt(spec)

# define usage function
print_usage <- function (spec = NULL) {
	cat(getopt(spec, usage = T))
	cat("Usage example: \n")
	cat("
Usage:
	--help 	NULL 	get this help
	--result 	character 	the input snp file [forced]
	--chr 	character 	the chr id to draw, or the chr name if been provided in the second column [forced]
	--threshold 	numeric 	the threshold
	--output 	character 	the output file prefix, default: pop
	--pcol 	character 	the col line for point
	--lcol 	character 	the col line for line
	--qcol 	character 	the pcol line for threshold
	--xlab 	character 	x axis title
	--ylab 	character 	y axis title, only can used: deltaa-index, loess, Gprime, ED, ridit [forced]
	--chrgap 	numeric 	gap size between chr in plot, default: 10% of average chr length
\n")
	q(status = 1)
}

##ED
#opt$result <- "pop.sliding.result";opt$chr <- "chr.list";opt$output <- "test.ED";opt$pcol <- "slidingED";opt$threshold <- 0.9995;opt$ylab <- "ED"

##index
#opt$result <- "pop.bootstrap.result";opt$chr <- "chr.list";opt$output <- "test.index";opt$pcol <- "delta";opt$lcol <- "slidingD";opt$qcol <- "CI";opt$ylab <- "delta-index"

##loess
#opt$result <- "pop.bootstrap.result";opt$pcol <- "delta";opt$qcol <- "CI";opt$lcol <- "loess";opt$ylab <- "loess";opt$output <- "test.loess";opt$chr <- "chr.list"

##Gprime
#opt$result <- "pop.Gprime";opt$pcol <- "G";opt$lcol <- "Gprime";opt$qcol <- "GprimeT";opt$ylab <- "Gprime";opt$output <- "test.Gprime";opt$chr <- "chr.list"

if (!is.null(opt$help)) {print_usage(spec)}
if (is.null(opt$result)) {print_usage(spec)}
if (is.null(opt$chr)) {print_usage(spec)}
if (is.null(opt$ylab)) {print_usage(spec)}
if (is.null(opt$output)) {opt$output <- "pop"}
if (!is.null(opt$threshold)) {opt$threshold <- as.numeric(opt$threshold)}
#if (is.null(opt$pcol)) {opt$pcol <- "delta"}
#if (is.null(opt$lcol)) {opt$lcol <- "slidingD"}
#if (is.null(opt$qcol)) {opt$qcol <- "CI"}
#if (is.null(opt$xlab)) {opt$xlab <- "chromosome"}
#if (is.null(opt$threshold)) {opt$threshold <- 0.9995}
#########################################################################################
chr <- read.table(opt$chr)
pos <- read.table(opt$result, head = T, comment.char = "^")
if(!("slidingI1" %in% colnames(opt$result))){opt$mutmap=T}
if ("pos1" %in% colnames(pos)) {
    pos$pos1 <- as.numeric(pos$pos1)
    pos$pos2 <- as.numeric(pos$pos2)
    pos$pos <- (pos$pos1 + pos$pos2) / 2
}

dfpos <- data.frame(CHR = pos$X.chr, BP = pos$pos, Delta = pos[[opt$pcol]])

if (!is.null(opt$lcol)) {
  dfpos$Slide <- pos[[opt$lcol]]
}

if (!is.null(opt$qcol)) {
  dfpos$CI <- pos[[opt$qcol]]
}

dfpos <- dfpos %>% filter(CHR %in% chr$V1)
#Gprime=pos$Gprime,G=pos$G,ED=pos$ED,EDprime=pos$EDprime,Gfdr=pos$Gfdr,EDfdr=pos$EDfdr)
if (ncol(chr) > 1) {
  dfpos$CHR <- chr[match(dfpos$CHR, chr$V1), 2]
  rownames(chr) <- chr$V2
} else {
  rownames(chr) <- chr$V1
}

lev <- NULL
lev$CHR <- levels(as.factor(dfpos$CHR))
lev$order <- gsub("\\D", "", lev$CHR)
lev$order <- as.numeric(lev$order)

dfpos <- merge(dfpos, lev, by = "CHR")
dfpos <- arrange(dfpos, order, BP)

if (is.null(opt$chrgap)) {
  opt$chrgap <- dfpos %>%
    group_by(order) %>%
    summarise(chr_len = max(BP)) %>%
    pull(chr_len) %>%
    mean() / 10}

dpos <- dfpos %>%
  group_by(order) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(as.numeric(chr_len) + as.numeric(opt$chrgap)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(dfpos, ., by = c("order" = "order")) %>%
  arrange(order, BP) %>%
  mutate(BPcum = BP + tot)

print("haha")

axisdf <- dpos %>%
  group_by(CHR) %>%
  summarize(center = (as.numeric(max(BPcum)) + as.numeric(min(BPcum))) / 2 )
#########################################################################################
# plot total chr
x_axis_start <- as.numeric(axisdf[which(axisdf$CHR == unique(dpos$CHR)[1]), 2])
x_axis_end <- as.numeric(axisdf[which(axisdf$CHR == unique(dpos$CHR)[length(unique(dpos$CHR))]), 2])
print(head(dpos))
plot_bsa <- ggplot(dpos) +
  geom_point(aes(x = BPcum, y = Delta, color = as.factor(order)), size = 1, alpha = 1) +
   scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
#  facet_grid(.~as.factor(order), switch = "x", space = "free", scales = "free_x") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_color_manual(values = (rep(c("slategray3", "skyblue4"), ceiling(nrow(chr) / 2)))[1:nrow(chr)]) +
  geom_segment(x = x_axis_start, xend = x_axis_end, y = -Inf, yend = -Inf) +
  coord_cartesian(clip = "off")

plot_bsa_ridit <- ggplot(dpos) +
  geom_point(aes(x = BPcum, y = Delta, color = as.factor(order))) +
   scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
#  facet_grid(.~as.factor(order), switch = "x", space = "free", scales = "free_x") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  scale_color_manual(values = (rep(c("slategray3", "skyblue4"), ceiling(nrow(chr) / 2)))[1:nrow(chr)]) +
  geom_segment(x = x_axis_start, xend = x_axis_end, y = -Inf, yend = -Inf) +
  coord_cartesian(clip = "off")

if (opt$ylab == "delta-index") {
#  plot_bsa <- plot_bsa + labs(y = expression(Delta*" (SNP-index)"))
#  plot_bsa <- plot_bsa + labs(y = TeX("$\\Delta (SNP-index)$"))
#  plot_bsa <- plot_bsa + labs(y = "Δ (SNP-index)")
  if (min(dpos$Delta) < 0) {
    plot_bsa <- plot_bsa + labs(y = "Δ SNP-index",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(-1, 1))
  } else if(opt$mutmap){
    plot_bsa <- plot_bsa + labs(y = "SNP-index",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(0, 1))
  }else{
    plot_bsa <- plot_bsa + labs(y = "|Δ SNP-index|",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(0, 1))
  }
} else if (opt$ylab == "loess") {
  if(opt$mutmap){
  plot_bsa <- plot_bsa + labs(y = "SNP-index",title="manhattan plot of index-loess for BSA")
  }else if(opt$abs){
    plot_bsa_chr <- plot_bsa_chr + labs(y = "|Δ SNP-index|",title="manhattan plot of index-loess for BSA")
  }else{
  plot_bsa <- plot_bsa + labs(y = "Δ SNP-index",title="manhattan plot of index-loess for BSA")
 }
  newdata<-na.omit(dpos$Delta)
  if (min(newdata) < 0) {
    plot_bsa <- plot_bsa + scale_y_continuous(limits = c(-1, 1))
  } else {
    plot_bsa <- plot_bsa + scale_y_continuous(limits = c(0, 1))
  }
} else if (opt$ylab == "ED") {
  plot_bsa <- plot_bsa + labs(y = "ED4",title="manhattan plot of ED for BSA") + scale_y_continuous(limits = c(0, max(dpos$Delta)*1.2))
} else if (opt$ylab == "Gprime") {
  plot_bsa <- plot_bsa + labs(y = "Gprime",title="manhattan plot of Gprime for BSA") + scale_y_continuous(limits = c(0, max(dpos$Delta)))
} else if (opt$ylab == "ridit") {
  plot_bsa <- plot_bsa_ridit + labs(y = '-log10(p)',title="manhattan plot of ridit for BSA")
  #plot_bsa <- plot_bsa + labs(y = expression(paste(italic(p)*" value"))) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
} else if (opt$ylab == "deepBSA_DL"){
  plot_bsa <- plot_bsa + labs(y = "DL_values",title="manhattan plot of DeepBSA_DL for BSA") + scale_y_continuous(limits = c(0, max(dpos$Delta)))
} else if (opt$ylab == "deepBSA_K"){
  plot_bsa <- plot_bsa + labs(y = "K_values",title="manhattan plot of DeepBSA_K for BSA") + scale_y_continuous(limits = c(0, max(dpos$Delta)))
}

if (!is.null(opt$lcol)) {
  plot_bsa <- plot_bsa + geom_line(mapping = aes(x = BPcum, y = Slide, group = order), color = "grey0")
}

if (!is.null(opt$qcol)) {
  plot_bsa <- plot_bsa + geom_line(mapping = aes(x = BPcum, y = CI, group = order), color = "red", linetype = "dashed")
#  if (opt$ylab == "delta-index") {
#    plot_bsa <- plot_bsa + geom_line(mapping = aes(x = BPcum, y = -1 * CI), color = "red")
#  }
}

if (!is.null(opt$threshold)) {
  quantile <- quantile(dfpos$Delta, opt$threshold, na.rm = T)
  plot_bsa <- plot_bsa + geom_hline(yintercept = quantile, linetype = "dashed", col = "red")
}

plot_bsa <- plot_bsa +
  theme(axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title=element_text(size = 12,face="bold",hjust=0.5),
#        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, margin = margin(0, 0.5, 0, 0, "cm")),
        axis.line.y  = element_line(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "none")

if (!is.null(opt$xlab)) {
  plot_bsa <- plot_bsa + labs(x = opt$xlab) + theme(axis.title.x = element_text(size = 15, margin = margin(0.5, 0, 0, 0, "cm")))
} else {
  plot_bsa <- plot_bsa + theme(axis.title.x = element_blank())
}

#ggsave(plot_bsa, file = paste(opt$out, ".index.pdf", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = "pdf")
ggsave(plot_bsa, file = paste(opt$out, ".index.pdf", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = cairo_pdf)
ggsave(plot_bsa, file = paste(opt$out, ".index.png", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = "png")
#########################################################################################
# plot single chr
#i <- rownames(chr)[2]
for (i in rownames(chr)) {
  plot_data_chr <- dpos[which(dpos$CHR == i), ]
#  plot_data_chr$BP <- as.numeric(plot_data_chr$BP)


  plot_bsa_chr <- ggplot(plot_data_chr) +
    geom_point(aes(x = BP, y = Delta), color = "slategray3", size = 1, alpha = 1) +
#    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
#    facet_grid(.~as.factor(order), switch = "x", space = "free", scales = "free_x") +
    guides(x = "axis_truncated", y = "axis_truncated") +
    geom_segment(x = min(plot_data_chr$BP), xend = max(plot_data_chr$BP), y = -Inf, yend = -Inf) +
    coord_cartesian(clip = "off")
  if (opt$ylab == "delta-index") {
#    plot_bsa_chr <- plot_bsa_chr + labs(y = expression(Delta*" (SNP-index)"), x = i)
#    plot_bsa_chr <- plot_bsa_chr + labs(y = TeX("$\\Delta (SNP-index)$"), x = i)
#    plot_bsa_chr <- plot_bsa_chr + labs(y = "Δ (SNP-index)", x = i)
     if (opt$abs) {
    plot_bsa_chr <- plot_bsa_chr + labs(y = "|Δ SNP-index|",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(0, 1))
  } else if(opt$mutmap){
    plot_bsa_chr <- plot_bsa_chr + labs(y = "SNP-index",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(0, 1))
  }else{
    plot_bsa_chr <- plot_bsa_chr + labs(y = "Δ SNP-index",title="manhattan plot of index-slid for BSA") + scale_y_continuous(limits = c(-1, 1))
  }
  } else if (opt$ylab == "loess") {
    if(opt$mutmap){
  plot_bsa_chr <- plot_bsa_chr + labs(y = "SNP-index",title="manhattan plot of index-loess for BSA")
    }else if(opt$abs){
  plot_bsa_chr <- plot_bsa_chr + labs(y = "|Δ SNP-index|",title="manhattan plot of index-loess for BSA")
 }else{
     plot_bsa <- plot_bsa_chr + labs(y = "Δ SNP-index",title="manhattan plot of index-loess for BSA")

 }
    newdata<-na.omit(dpos$Delta)
    if (min(newdata) < 0) {
      plot_bsa_chr <- plot_bsa_chr + scale_y_continuous(limits = c(-1, 1))
    } else {
      plot_bsa_chr <- plot_bsa_chr + scale_y_continuous(limits = c(0, 1))
    }
  } else if (opt$ylab == "ED") {
    plot_bsa_chr <- plot_bsa_chr + labs(y = "ED", x = i) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
  } else if (opt$ylab == "Gprime") {
    plot_bsa_chr <- plot_bsa_chr + labs(y = "Gprime", x = i) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
  } else if (opt$ylab == "ridit") {
    plot_bsa <- plot_bsa_ridit + labs(y = '-log10(p)',title="manhattan plot of ridit for BSA")
    #plot_bsa_chr <- plot_bsa_chr + labs(y = expression(paste(italic(p)*" value")), x = i) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
  } else if (opt$ylab == "deepBSA_DL"){
    plot_bsa_chr <- plot_bsa_chr + labs(y = "DL_values", x = i) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
  } else if (opt$ylab == "deepBSA_K"){
    plot_bsa_chr <- plot_bsa_chr + labs(y = "K_values", x = i) + scale_y_continuous(limits = c(0, max(dpos$Delta)))
  }

  if (!is.null(opt$lcol)) {
    plot_bsa_chr <- plot_bsa_chr + geom_line(mapping = aes(x = BP, y = Slide), color = "grey0")
  }

  if (!is.null(opt$qcol)) {
    plot_bsa_chr <- plot_bsa_chr + geom_line(mapping = aes(x = BP, y = CI), color = "red", linetype = "dashed")
#    if (opt$ylab == "delta-index") {
#      plot_bsa_chr <- plot_bsa_chr + geom_line(mapping = aes(x = BP, y = -1 * CI), color = "red")
#    }
  }
# write.table(dfpos,"dfpos.txt")
  if (!is.null(opt$threshold)) {
    quantile <- quantile(dfpos$Delta, opt$threshold, na.rm = T)
    plot_bsa_chr <- plot_bsa_chr + geom_hline(yintercept = quantile, linetype = "dashed", col = "red")
  }

  plot_bsa_chr <- plot_bsa_chr +
    theme(axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text = element_text(size = 12, color = "black"),
#          axis.text.x = element_text(angle = 45, hjust = 1),
#          axis.title.x = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(0.5, 0, 0, 0, "cm")),
          axis.title.y = element_text(size = 15, margin = margin(0, 0.5, 0, 0, "cm")),
#          axis.line.y  = element_line(),
          axis.line  = element_line(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0, "lines"),
          legend.position = "none")

#  ggsave(plot_bsa_chr, file = paste(opt$out, ".", chr[i, 1], ".index.pdf", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = "pdf")
  ggsave(plot_bsa_chr, file = paste(opt$out, ".", chr[i, 1], ".index.pdf", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = cairo_pdf)
  ggsave(plot_bsa_chr, file = paste(opt$out, ".", chr[i, 1], ".index.png", sep = ""), width = 14, height = 6, dpi = 300, units = "in", device = "png")
}
#########################################################################################
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
