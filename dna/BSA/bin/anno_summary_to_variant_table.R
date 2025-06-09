#!/usr/bin/env Rscript
# @Last-edit Time 2023/9/26
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com
options(bitmapType = "cairo")
library(getopt)
library(data.table)
library(stringr)
library(dplyr)
spec <- matrix(c(
    "vt", "v", 1, "character",
    "anno", "a", 2, "character",
    "outname", "o", 1, "character",
    "wgs", "w", 2, "character",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec = spec) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("
Usage:
        --vt,-v   variant table
        --anno,-a   anno summary, 提供wgs时可不填
        --outname,-o   输出文件名
        --wgs,-w   wgsv4工作目录，用来抓取anno summary，可不填
        --help,-h   usage
\n")
    q(status = 1)
}
if (!is.null(opt$help)) {
    print_usage(spec)
}
if (is.null(opt$vt)) {
    print_usage(spec)
}
if (is.null(opt$outname)) {
    print_usage(spec)
}
if (!is.null(opt$anno)) {
    anno <- opt$anno
} else if (!is.null(opt$wgs)) {
    anno <- file.path(opt$wgs, "output/tmp/02.reference/anno.summary")
} else {
    print_usage(spec)
}
anno_df <- as.data.frame(fread(anno, sep = "\t", header = TRUE))
id_df <- as.data.frame(str_match(anno_df$GeneID, "(\\S+):(\\S+):(\\S+):(\\S+):(\\S+)"))
names(id_df) <- c("GeneID", "Gene_id", "Transcript_id", "Chr", "Pos1", "Pos2")
anno_df <- cbind(anno_df[, names(anno_df) != "GeneID"], id_df[, -1])
vt_df <- as.data.frame(fread(opt$vt, sep = "\t", header = TRUE))
vt_df$Transcript_id <- vt_df$Feature_id
feature_in <- vt_df$Transcript_id %in% anno_df$Transcript_id
feature_abs_in <- gsub("\\.\\d+$", "", vt_df$Transcript_id) %in% anno_df$Transcript_id
vt_df$Transcript_id[(!feature_in) & feature_abs_in] <- gsub("\\.\\d+$", "", vt_df$Transcript_id[(!feature_in) & feature_abs_in])
vt_df <- vt_df %>% left_join(anno_df[, !names(anno_df) %in% c("Gene_id", "Chr", "Pos1", "Pos2")], by = "Transcript_id")
## replacement <- rep("--", ncol(vt_df))
## names(replacement) <- names(vt_df)
## vt_df <- tidyr::replace_na(vt_df, as.list(replacement))
write.table(vt_df, opt$outname, sep = "\t", quote = FALSE, row.names = FALSE, na = "--")
print("Done!")
