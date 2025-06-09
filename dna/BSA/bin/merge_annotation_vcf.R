#!/usr/bin/env Rscript

if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)
###传参信息
spec <- matrix(c(
  'variant_table', 'r', 0, 'character',
  'anno_summary', 'l', 0, 'character',
  'outfile', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
Usage:
	--variant_table    输入variant_table文件
	--anno_summary     输入anno_summary列表
	--outfile   输出文件
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$variant_table))   { print_usage(spec)}
if ( is.null(opt$anno_summary))  { print_usage(spec) }
if ( is.null(opt$outfile))  { print_usage(spec) }

pacman::p_load(tidyverse, purrr, data.table)

### 导入variant.table文件
variant_table <-  "/mnt/ilustre/users/yuan.xu/majorbio_project/MJ20220708119-MJ-R-20220719043-燕景慧-WGS-BSA-4个样本-20230217-FX2023021700141/workflow_results/04.snpIndel/finalvcf/variant.table"
variant_table <- fread(opt$variant_table, header = TRUE, sep="\t", fill=TRUE)

### 导入anno_summary文件
anno_summary <- "/mnt/ilustre/users/yuan.xu/majorbio_project/MJ20220708119-MJ-R-20220719043-燕景慧-WGS-BSA-4个样本-20230217-FX2023021700141/workflow_results/02.reference/anno.summary"
anno_summary <- fread(opt$anno_summary, header = TRUE, sep="\t", fill=TRUE)
new_anno_summary <- anno_summary %>%
  separate(col = GeneID, into = c("Gene_Name", "Gene_ID", "CHROM", "gene_start", "gene_end"), sep = ":")

result_df <- new_anno_summary %>%
  right_join(variant_table, by.y = "Feature_ID", by.x = "Gene_ID") %>%
  select(CHROM, POS, REF, Gene_Name, Gene_ID, gene_start, gene_end, NRID, NRANNO, UniID, UniANNO, KoID, Koanno, GOTERM,
         GOANNO, EGGNOG, EGGNOG_ANNO, PfamID, PfamAnno, InterProAccession, InterProAnno, ANN_Annotation, Annotation_impact,
         Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HVGS.p)

###结果输出
write.table(result_df, opt$outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)