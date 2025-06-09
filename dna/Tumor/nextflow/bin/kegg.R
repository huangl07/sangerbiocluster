#!/usr/bin/env Rscript
times<-Sys.time()

if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(getopt)

###传参信息
spec <- matrix(c(
  'degfile', 'd', 0, 'character',
  'term2genefile', 'g', 0, 'character',
  'term2namefile', 'n', 0, 'character',
  'outname', 'o', 0, 'character',
  'db','b','0','character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--degfile    输入基因列表列表
	--term2genefile   输入term2gene文件
  --term2namefile   输入term2name文件
	--outname   输出文件名字
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$degfile))   { print_usage(spec)}
if ( is.null(opt$term2genefile))  { print_usage(spec) }
if ( is.null(opt$term2namefile))  { print_usage(spec) }
if ( is.null(opt$outname))  { print_usage(spec) }

library(org.Ddemo.eg.db,lib.loc=opt$db)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)

eg <- read.table(opt$degfile,sep = '\t',header = F, quote = '')[,1]
DEG <- bitr(eg,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb='org.Hs.eg.db')
term2gene <- read.table(opt$term2genefile,sep = '\t',header = F, quote = '')
term2name <- read.table(opt$term2namefile,sep = '\t',header = F, quote = '')
Species_db <- "org.Ddemo.eg.db"
outname <- opt$outname

###GO富集
dir.create('GO_result', showWarnings = FALSE)
GO_enrich <- enrichGO(    gene         = DEG$ENTREZID,
                          OrgDb         = "org.Hs.eg.db",
                          keyType       = 'ENTREZID',
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 1,
                          qvalueCutoff  = 1,
                          readable      = FALSE,minGSSize=0)
GO_enrich_df <- as.data.frame(GO_enrich)
GO_enrich_df <- GO_enrich_df[order(GO_enrich_df$p.adjust),] #按照p.adjust排序

write.table(GO_enrich_df,paste0('GO_result/',outname, "_GOenrichment.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)

#做富集气泡图
GO_enrich_df <- GO_enrich_df[1:20,]
GO_enrich_df <- GO_enrich_df[order(GO_enrich_df$ONTOLOGY),]
GO_enrich_df$Description <- factor(GO_enrich_df$Description, levels = GO_enrich_df$Description )
GO_enrich_df <- GO_enrich_df %>% 
  mutate(GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio)))
p1 <- ggplot(GO_enrich_df, aes(x=GeneRatio_result,y=reorder(Description, GeneRatio_result))) +
  geom_point(aes(size=Count, color=p.adjust, shape=ONTOLOGY)) +
  coord_cartesian(clip="off") + 
  scale_color_gradient(low="red", high="blue",name="p.adjust") +
  theme_bw() +
  scale_size_continuous(name="Count",range=c(1,10))+
  labs(x="GeneRatio", y="Terms")
ggsave(p1, file = file.path(paste0('GO_result/',outname, "_GOenrichment.png")), dpi = 600, height = 10, width = 15)
ggsave(p1, file = file.path(paste0('GO_result/',outname, "_GOenrichment.pdf")), dpi = 600, height = 10, width = 15)

###KEGG富集

dir.create('KEGG_result', showWarnings = FALSE)
KEGG_enrich <- enricher(gene = DEG$SYMBOL,
                      TERM2GENE = term2gene,
                      TERM2NAME = term2name,
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,minGSSize=0)
KEGG_enrich_0.05 <- enricher(gene = DEG$SYMBOL,
                      TERM2GENE = term2gene,
                      TERM2NAME = term2name,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 1,minGSSize=0)
KEGG_enrich_df <- as.data.frame(KEGG_enrich)
KEGG_enrich_df
KEGG_enrich_df <- KEGG_enrich_df[order(KEGG_enrich_df$p.adjust),] #按照p.adjust排序
KEGG_enrich_0.05_df <- as.data.frame(KEGG_enrich_0.05)
KEGG_enrich_0.05_df <- KEGG_enrich_0.05_df[order(KEGG_enrich_0.05_df$p.adjust),] #按照p.adjust排序
write.table(KEGG_enrich_df,paste0('KEGG_result/',outname, "_KEGGenrichment.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
write.table(KEGG_enrich_0.05_df,paste0('KEGG_result/',outname, "_KEGGenrichment_0.05.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
#做富集气泡图
KEGG_enrich_df <- KEGG_enrich_df[1:20,]
KEGG_enrich_df$Description <- factor(KEGG_enrich_df$Description, levels = KEGG_enrich_df$Description )
KEGG_enrich_df <- KEGG_enrich_df %>% 
mutate(GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio)))
p2 <- ggplot(KEGG_enrich_df, aes(x=GeneRatio_result,y=reorder(Description, GeneRatio_result))) +
geom_point(aes(size=Count, color=p.adjust)) +
coord_cartesian(clip="off") + 
scale_color_gradient(low="red", high="blue",name="p.adjust") +
theme_bw() +
scale_size_continuous(name="Count",range=c(1,10))+
labs(x="GeneRatio", y="Terms")
ggsave(p2, file = file.path(paste0('KEGG_result/',outname, "_KEGGenrichment.png")), dpi = 600, height = 10, width = 15)
ggsave(p2, file = file.path(paste0('KEGG_result/',outname, "_KEGGenrichment.pdf")), dpi = 600, height = 10, width = 15)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
