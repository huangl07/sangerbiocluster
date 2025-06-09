#!/usr/bin/env Rscript
times <- Sys.time()
library(getopt)
###传参信息
spec <- matrix(c(
  'dbnfspfile', 'd', 0, 'character',
  'mergeannofile', 'm', 0, 'character',
  'outmergefile', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--dbnfspfile        输入dbNFSP的表格
	--mergeannofile     输入merge_anno的表格
	--outmergefile      输出的肿瘤易感基因的详细信息
	--help      usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$dbnfspfile))   { print_usage(spec)}
if ( is.null(opt$mergeannofile))  { print_usage(spec) }
if ( is.null(opt$outmergefile))  { print_usage(spec) }

library(tidyverse)

dbNFSP_result <- read.delim(opt$dbnfspfile, header = FALSE, sep = "\t")
colnames(dbNFSP_result) <- c('chromosome', 'position', 'rs_id', 'ref', 'alt', 'MetaSVM', 'FATHMM', 'LRT', 'PROVEAN', 'MutationTaster', 'MutationAssessor', 'Polyphen2_HDIV', 'Polyphen2_HVAR', 'SIFT')

merge_anno_result <- read.delim(opt$mergeannofile, header = TRUE, sep = "\t")

prediction_info <- dbNFSP_result %>% 
  mutate(MetaSVM_pred = map_chr(strsplit(MetaSVM, ","), ~ case_when("D" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(FATHMM_pred = map_chr(strsplit(FATHMM, ","), ~ case_when("D" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(LRT_pred = map_chr(strsplit(LRT, ","), ~ case_when("D" %in% .x ~ "D", TRUE ~ "T"))) %>%
  mutate(PROVEAN_pred = map_chr(strsplit(PROVEAN, ","), ~ case_when("D" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(MutationTaster_pred = map_chr(strsplit(MutationTaster, ","), ~ case_when("A" %in% .x ~ "D", "D" %in% .x ~ "D",TRUE ~ "T"))) %>% 
  mutate(MutationAssessor_pred = map_chr(strsplit(MutationAssessor, ","), ~ case_when("H" %in% .x ~ "D", "M" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(Polyphen2_HDIV_pred = map_chr(strsplit(Polyphen2_HDIV, ","), ~ case_when("D" %in% .x ~ "D", "P" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(Polyphen2_HVAR_pred = map_chr(strsplit(Polyphen2_HVAR, ","), ~ case_when("D" %in% .x ~ "D", "P" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  mutate(SIFT_pred = map_chr(strsplit(SIFT, ","), ~ case_when("D" %in% .x ~ "D", TRUE ~ "T"))) %>% 
  select(chromosome, position, rs_id, ref, alt, MetaSVM_pred, FATHMM_pred, LRT_pred, PROVEAN_pred, MutationTaster_pred, MutationAssessor_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, SIFT_pred) %>% 
  mutate(Deleterious_count = rowSums(. == "D"))
  
Deleterious_site <- prediction_info %>% 
  filter(Deleterious_count >= 3) %>% 
  merge(merge_anno_result, by = c('chromosome', 'position', 'ref', 'alt')) %>% 
  select(-dbsnp_id)

write.table(Deleterious_site, file = opt$outmergefile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

escape_time <- Sys.time()-times
print("Done!")
