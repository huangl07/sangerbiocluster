#!/usr/bin/env Rscript
times <- Sys.time()
library("getopt")
options(bitmapType = "cairo")
options(scipen = 200)
spec <- matrix(c(
    "infile", "i", 0, "character",
    "outfile1", "o", 0, "character",
    "outfile2", "v", 0, "character",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec = NULL) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("
Usage:
	--infile	the input file
	--outfile1	the output file(png)
    --outfile2  the output file(pdf)
	--help		usage
\n")
    q(status = 1)
}
if (!is.null(opt$help)) {
    print_usage(spec)
}
if (is.null(opt$infile)) {
    print_usage(spec)
}
if (is.null(opt$outfile1)) {
    print_usage(spec)
}
if (is.null(opt$outfile2)) {
    print_usage(spec)
}
pacman::p_load(data.table, ggseqlogo, tidyverse, patchwork)

grep_string <- paste0("grep NUCLEOTIDE ", opt$infile)
# grep_string = paste0("grep NUCLEOTIDE ", "/mnt/lustre/users/sanger-dev/sg-users/xuyuan/majorbio_development/crispr_off_target/20230818_test/demo/01.homo_region/homo_region.extended_profile.xls")
profile <- as.data.frame(fread(cmd = grep_string))
profile <- profile[, -ncol(profile)]
A_sum <- apply(profile[profile$V1 == "NUCLEOTIDE A", -1], 2, sum)
C_sum <- apply(profile[profile$V1 == "NUCLEOTIDE C", -1], 2, sum)
G_sum <- apply(profile[profile$V1 == "NUCLEOTIDE G", -1], 2, sum)
T_sum <- apply(profile[profile$V1 == "NUCLEOTIDE T", -1], 2, sum)
sum_data <- rbind(A_sum, C_sum, G_sum, T_sum)
rownames(sum_data) <- c("A", "C", "G", "T")
p1 <- ggseqlogo(sum_data, method = "bits")
p2 <- ggseqlogo(sum_data, method = "prob")
g <- p1 / p2
ggsave(
    filename = opt$outfile1, # 保存的文件名称。通过后缀来决定生成什么格式的图片
    plot = g,
    width = 16, # 宽
    height = 9, # 高
    units = "in", # 单位
    dpi = 300, # 分辨率DPI
    bg = "white"
)
ggsave(
    filename = opt$outfile2, # 保存的文件名称。通过后缀来决定生成什么格式的图片
    plot = g,
    width = 16, # 宽
    height = 9, # 高
    units = "in", # 单位
    dpi = 300, # 分辨率DPI
    bg = "white"
)
