#!/usr/bin/env Rscript
###times <- Sys.time()
###library(getopt)
######传参信息
###spec <- matrix(c(
###  'seqgz', 's', 0, 'character',
###  'outpath', 'o', 0, 'character',
###  'samplename', 'n', 0, 'character',
###  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
###opt <- getopt(spec)
###print_usage <- function(spec=NULL){
###  cat(getopt(spec, usage=TRUE));
###  cat("Usage example: \n")
###  cat("	
###Usage:
###	--seqgz    输入data.bin50_seqz.gz
###	--outpath
###    --samplename
###	--help    usage
###\n")
###  q(status = 1);
###}
###if ( !is.null(opt$help))   { print_usage(spec) }
###if ( is.null(opt$seqgz))   { print_usage(spec)}
###if ( is.null(opt$outpath))  { print_usage(spec) }
###if ( is.null(opt$samplename))  { print_usage(spec) }

library(GenomicRanges)
library(magrittr)
#准备实例数据
snps <- read.delim("./Test_combine_2.txt", head=TRUE, sep="\t")
cnvs <- read.delim("./Test_combine_1.txt", head=TRUE, sep="\t")
gsnps <- GRanges(seqnames = snps$chr ,
               ranges = IRanges(snps$start , snps$end ),
               strand = "+" )
#metadata columns can be added to a GRanges object,将表格中的其他信息添加到gsnps中
mcols(gsnps) <- snps

gcnvs <- GRanges(seqnames = cnvs$chr,
               ranges = IRanges(cnvs$start , cnvs$end ),
               strand = "+" )
#metadata columns can be added to a GRanges object
mcols(gcnvs) <- cnvs
overlaps <- findOverlaps(gsnps, gcnvs)
overlaps

#合并信息。
merge <- cbind( mcols(gsnps[queryHits(overlaps), ]) , mcols(gcnvs[subjectHits(overlaps) ,]) )

#整合信息
merge <- as.data.frame(merge) %>%
     dplyr::mutate(mutation_id = rs,
                normal_cn = 2,
                ) %>%
  dplyr::select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)

write.table(merge, file ="pyclone_input.tsv", sep="\t", row.names =TRUE, 
            col.names =FALSE, quote =FALSE)