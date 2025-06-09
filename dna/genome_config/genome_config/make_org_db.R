# @Last-edit Time 2022/12/1
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com

library(tidyverse)
library(AnnotationForge)
library(getopt)

spec <- matrix(c(
    "annofile", "a", 1, "character",
    "outprefix", "o", 1, "character",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print_usage <- function(spec = NULL) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("
Usage:
    --annofile,-a     输入整理好的GO和KEGG注释列表
    --outprefix,-o    输出orgdb构建路径
    --help,-h         usage
\n")
    q(status = 1)
}
if (!is.null(opt$help)) {
    print_usage(spec)
}
if (is.null(opt$annofile)) {
    print_usage(spec)
}
if (is.null(opt$outprefix)) {
    print_usage(spec)
}

anno_result <- read_delim(opt$annofile, delim = "\t", quote = "")
gene_name <- anno_result %>%
    dplyr::select(transcript_id) %>%
    distinct()

gene_info <- tibble(
    GID = gene_name$transcript_id,
    Gene_Name = gene_name$transcript_id
)

gene2go <- anno_result %>%
    dplyr::select(transcript_id, go_term) %>%
    as_tibble() %>%
    filter(go_term != "", go_term != "--") %>%
    separate_rows(go_term, sep = ",", convert = FALSE) %>%
    distinct(transcript_id, go_term, .keep_all = TRUE) %>%
    filter(grepl("^GO:", go_term)) %>%
    mutate(EVIDENCE = "IEA")
colnames(gene2go) <- c("GID", "GO", "EVIDENCE")

f_chr <- tibble(
    GID = anno_result$transcript_id,
    CHROMOSOME = anno_result$chr
) %>%
    distinct(GID, .keep_all = TRUE)

if (!dir.exists(opt$outprefix)) {
    dir.create(opt$outprefix)
}

makeOrgPackage(
    gene_info = gene_info,
    chromosome = f_chr,
    go = gene2go,
    maintainer = "XY<yuan.xu@majorbio.com>",
    author = "yuan.xu",
    version = "1.0",
    outputDir = opt$outprefix,
    tax_id = "1",
    genus = "demo",
    species = "demo",
    goTable = "go"
)

term2gene <- anno_result %>%
    dplyr::select(ko_id, transcript_id) %>%
    filter(ko_id != "--", ko_id != "") %>%
    distinct() %>%
    separate_rows(ko_id, sep = ",", convert = FALSE) %>%
    group_by(ko_id) %>%
    nest()

gsid2gene <- lapply(seq_len(nrow(term2gene)), function(i) {
    term2gene$data[[i]]$transcript_id
})
names(gsid2gene) <- term2gene$ko_id

gsid2name <- anno_result %>%
    dplyr::select(ko_id, ko_anno) %>%
    filter(ko_id != "--", ko_id != "") %>%
    mutate(ko_id = gsub(pattern = ",", replacement = ":", .$ko_id, )) %>%
    separate_rows(ko_id, ko_anno, sep = ":") %>%
    distinct(ko_id, .keep_all = TRUE) %>%
    as.list()
names(gsid2name) <- c("gsid", "name")

gson <- list(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    gene2name = list(),
    species = c("demo"),
    gsname = c("KEGG"),
    version = c("1.0"),
    accessed_date = c(format(Sys.Date(), "%Y-%m-%d")),
    info = c()
)
jsonlite::write_json(gson, file.path(opt$outprefix, "kegg.gson"), pretty = TRUE)
