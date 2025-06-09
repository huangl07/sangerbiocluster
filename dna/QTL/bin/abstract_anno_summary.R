#!/usr/bin/env Rscript
times <- Sys.time()
if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)
###传参信息
spec <- matrix(c(
  'regionfile', 'r', 0, 'character',
  'anno_summary_file', 'g', 0, 'character',
  'variant_table_file', 'v', 0, 'character',
  'outfile1', 'o', 0, 'character',
  'outfile2', 'p', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--regionfile    输入pop.ED.region列表
	--anno_summary_file   输入anno_summary文件
  --variant_table_file variant文件
	--outfile1   输出文件1
  --outfile2   输出文件2
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$regionfile))   { print_usage(spec)}
if ( is.null(opt$anno_summary_file))  { print_usage(spec) }
if ( is.null(opt$variant_table_file))  { print_usage(spec) }
if ( is.null(opt$outfile1))  { print_usage(spec) }
if ( is.null(opt$outfile2))  { print_usage(spec) }

pacman::p_load(tidyverse, purrr)

 pop_region=read.table(opt$regionfile)
 if(grepl("_",pop_region[1,1])){
     pop_region <- pop_region %>%
   separate(V1, sep="[-_]", into = c("chr", "pos_start1", "pos_end1")) %>%
   separate(V2, sep="[-_]", into = c("chr", "pos_start2", "pos_end2"))
     pop_region$pos1 <- with(pop_region, pmin(as.numeric(pos_start1), as.numeric(pos_end1), as.numeric(pos_start2), as.numeric(pos_end2)))
     pop_region$pos2 <- with(pop_region, pmax(as.numeric(pos_start1), as.numeric(pos_end1), as.numeric(pos_start2), as.numeric(pos_end2)))
 }else{
     pop_region <- pop_region %>%
         separate(V1, sep="[-_]", into = c("chr", "pos_start1")) %>%
         separate(V2, sep="[-_]", into = c("chr", "pos_start2"))
     pop_region$pos1=with(pop_region,pmin(as.numeric(pos_start1),as.numeric(pos_start2)));
     pop_region$pos2=with(pop_region,pmax(as.numeric(pos_start1),as.numeric(pos_start2)))
 }
pop_region <- pop_region %>%
  select(chr, pos1, pos2)
print(pop_region)
### 导入anno_summary文件
#anno_summary_file <- "/mnt/lustre/users/sanger-dev/sg-users/yuan.xu/majorbio_task/遗传图谱报告测试/output/tmp/02.reference/anno.summary"
anno_summary <- read.delim(opt$anno_summary_file, header = TRUE) %>%
  as_tibble()
anno_summary_df <- anno_summary %>% 
  separate(GeneID, into=c("Gene_id", "Transcript_id", "Chrom", "transcript_start", "transcript_end"), sep=":", remove = FALSE)


###导入所有的gene_list

abstract_df <- function(x, pos_start, pos_end, anno_summary_df){
  df_abstract <- anno_summary_df %>%
    filter(Chrom == x, as.numeric(transcript_end) > as.numeric(pos_start), as.numeric(transcript_start) < as.numeric(pos_end)) %>%
    mutate(REGION = paste(x, pos_start, pos_end, sep=":")) %>%
    select(REGION, Chrom, transcript_start, transcript_end, Gene_id, Transcript_id, NRID, NRANNO,
    UniID, UniANNO, KoID, Koanno, GOTERM, GOANNO, EGGNOG, EGGNOG_ANNO, PfamID, PfamAnno)
  return(df_abstract)
}
df_abstract_all <- pmap_dfr(pop_region[1:3], ~abstract_df(..1, ..2, ..3, anno_summary_df))

### 导入variant_table

#variant_table_file <- "/mnt/lustre/users/sanger-dev/sg-users/yuan.xu/majorbio_task/遗传图谱报告测试/qtl_workdir6/variant.diff.table"
variant_table <- read.delim(opt$variant_table_file, header = TRUE, check.names = FALSE) %>%
  as_tibble()

abstract_variant_df <- function(x, pos_start, pos_end, variant_table){
  df_abstract <- variant_table %>%
    filter(CHROM == x, as.numeric(pos_start) <= as.numeric(POS), as.numeric(POS) <= as.numeric(pos_end))
  return(df_abstract)
}

df_abstract_variant_all <- pmap_dfr(pop_region[1:3], ~abstract_variant_df(..1, ..2, ..3, variant_table))

merge_df <- merge(df_abstract_variant_all, df_abstract_all, by.x="Feature_ID", by.y="Transcript_id") %>%
  as_tibble() %>%
  select(-Chrom)

###结果输出
write.table(df_abstract_all, opt$outfile1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(merge_df, opt$outfile2, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
