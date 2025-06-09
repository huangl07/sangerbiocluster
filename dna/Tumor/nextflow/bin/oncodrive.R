#!/usr/bin/env Rscript
start_time <- Sys.time()
library(getopt)
spec <- matrix(c(
  'maf', 'm', 0, 'character',
  'out', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--maf    输入maf表格
  --out    输出基因突变全景图
	--help   usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$maf))   { print_usage(spec) }
if ( is.null(opt$out))   { print_usage(spec) }
library(maftools)
library(ggrepel)
library(ggplot2)
options(bitmapType = "cairo")
vc_cols =  c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", 
             "#8491B4B2", "#91D1C2B2", "#DC0000B2")
names(vc_cols)=c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
laml<-read.maf(opt$maf)
maf<-laml@data
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
res = laml.sig
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < 0.05, yes = 'sig', no = 'nonsig')
res[, log_fdr := -log10(as.numeric(fdr))]
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
res$color = ifelse(test = res$fdr < 0.05, yes = colCode[1], no = colCode[2])
par(mar = c(4, 4, 2, 2))
p1<-ggplot(res)+ geom_point(aes(x=fract_muts_in_clusters,y=log_fdr, size=log_fdr, color=color))+ 
geom_text_repel(aes(x=fract_muts_in_clusters,y=log_fdr,label=label),max.overlaps =30,nudge_x=0.05)+
theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
theme_classic(base_size = 16)+labs(x='Fraction of variants within clusters', y='-log10(fdr)')
ggsave(paste0(opt$out,"_plot_oncodriven.png"),p1)
ggsave(paste0(opt$out,"_plot_oncodriven.pdf"),p1)
write.table(res, paste0(opt$out, "_oncogene.xls"), sep = "\t", quote = FALSE, col.names = NA)
