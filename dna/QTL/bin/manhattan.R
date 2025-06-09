times <- Sys.time()
options(bitmapType = "cairo")

library(getopt)
command <- matrix(
    c(
        "help", "h", 0, "loical", "显示此帮助文件",
        "input", "i", 1, "character", "输入文件",
        "output", "o", 2, "character", "输出文件"
    ),
    byrow = T, ncol = 5
)
args <- getopt(command)

if (!is.null(args$help) || is.null(args$input)) {
    cat(paste(getopt(command, usage = T), "\n"))
    q(status = 1)
}

if (is.null(args$output)) {
    args$output <- "./manha.pdf"
}
data <- read.table(args$input, head = T)
threshold <- data[2, 5]
Rdata <- data[, c(2:4)]

library(dplyr)
library(ggplot2)
library(stringr)
Rdata <- Rdata %>% arrange(str_rank(chr, numeric = TRUE), pos)
Rdata$chr <- factor(Rdata$chr,
    levels = str_sort(unique(Rdata$chr), numeric = TRUE)
)

# 计算chr长度
chr_len <- Rdata %>%
    group_by(chr) %>%
    summarise(chr_len = max(pos))


# 计算每条chr的初始位置
chr_pos <- chr_len %>%
    mutate(chr.start = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len)


# 计算累计SNP的位置
Snp_pos <- chr_pos %>%
    left_join(Rdata, ., by = "chr") %>%
    arrange(chr, pos) %>%
    mutate(Chromosome = pos + chr.start)


# 在每条chr的中间准备X轴标签位置
X_axis <- Snp_pos %>%
    group_by(chr) %>%
    summarize(center = (max(Chromosome) + min(Chromosome)) / 2)

print(Snp_pos)
p <- ggplot(Snp_pos, aes(x = Chromosome, y = lod)) +
    # 设置线
    geom_line(aes(color = as.factor(chr)), alpha = 0.8, size = 1) +
    # 设置颜色
    scale_color_manual(values = rep(c("slategray3", "skyblue4"), 30)) +
    # 设定X轴
    scale_x_continuous(label = X_axis$chr, breaks = X_axis$center) +
    # 设置y范围
    # scale_y_continuous(limits = c(0, 14)) +
    # ylim(0,14) +
    ylim(0, (max(Snp_pos$lod) / 0.8)) +
    # 删除图例
    # guides(as.factor(chr) = FALSE)
    # guides(colour=guide_legend(title = NULL)) +
    # 添加阈值线
    # geom_hline(yintercept = c(3, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'),  linetype = c("dotted", "twodash")) +
    geom_hline(yintercept = c(threshold), color = c("red"), linetype = c("twodash")) +

    # 设置主题
    theme_bw() +
    theme(
        legend.position = "none",
        # panel.border = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # legend.title=element_blank()
    )
out <- getwd()
ggsave(p, file = paste(args$out, "pdf", sep = "."), device = "pdf", width = 16, height = 9)
ggsave(p, file = paste(args$out, "png", sep = "."), device = "png", width = 16, height = 9)


escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
q()
