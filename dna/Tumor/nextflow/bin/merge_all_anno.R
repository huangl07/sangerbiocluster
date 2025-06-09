#!/usr/bin/env Rscript
times <- Sys.time()
library(getopt)
###传参信息
spec <- matrix(c(
  'cosmicfile', 'c', 0, 'character',
  'dbsnpfile', 'd', 0, 'character',
  'gnomADfile', 'g', 0, 'character',
  'snpefffile', 's', 0, 'character',
  'outmergefile', 'o', 0, 'character',
  'outvcffile', 'v', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--cosmicfile    输入cosmic的表格
	--dbsnpfile     输入dbsnp的表格
	--gnomADfile    输入gnomAD的表格
	--snpefffile    输入snpeff的表格
	--outmergefile  输出的联合分析结果
	--outvcffile    输出的vcf转换格式结果
	--help      usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$cosmicfile))   { print_usage(spec)}
if ( is.null(opt$dbsnpfile))  { print_usage(spec) }
if ( is.null(opt$gnomADfile))  { print_usage(spec) }
if ( is.null(opt$snpefffile))   { print_usage(spec)}
if ( is.null(opt$outmergefile))  { print_usage(spec) }
if ( is.null(opt$outvcffile))  { print_usage(spec) }

library(tidyverse)

cosmic_result <- read.delim(opt$cosmicfile, header = FALSE, sep = "\t")
colnames(cosmic_result) <- c('chromosome', 'position', 'cosmic_id', 'ref', 'alt', 'protein_change', 'aa_change')

dbsnp_result <- read.delim(opt$dbsnpfile, header = FALSE, sep = "\t")
colnames(dbsnp_result) <- c('chromosome', 'position', 'dbsnp_id', 'ref', 'alt', 'common') ##第6个字段common，若为1，则为common-snp

gnomAD_result <- read.delim(opt$gnomADfile, header = FALSE, sep = "\t")
colnames(gnomAD_result) <- c('chromosome', 'position', 'rs_id', 'ref', 'alt', 'AF_eas')

snpeff_result <- read.delim(opt$snpefffile, header = TRUE, sep = "\t")
colnames(snpeff_result) <- c('chromosome', 'position', 'ref', 'alt', 'effect', 'impact', 'gene', 'geneid')

###结果联合
merge_result <- merge(dbsnp_result, cosmic_result, by = c("chromosome", "position", "ref", "alt"), all = TRUE) %>% 
  filter(common != "1") %>% 
  filter(cosmic_id != "NA") %>% 
  merge(gnomAD_result, by = c("chromosome", "position", "ref", "alt")) %>% 
  select(-common, -rs_id) %>% 
  merge(snpeff_result, c("chromosome", "position", "ref", "alt"))
write.table(merge_result, file = opt$outmergefile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


merge_result2vcf <- merge_result %>% 
  select(`#CHROM` = "chromosome", POS = "position", ID = "dbsnp_id", REF = "ref", ALT = "alt") %>% 
  mutate(QUAL = ".", FILTER = "PASS", INFO = ".")
write.table(merge_result2vcf, file = opt$outvcffile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

