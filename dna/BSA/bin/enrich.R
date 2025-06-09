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
  'outdir','p','0','character',
  'title','t','0','character',
  'genename','a','0','character',
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
  --outdir  输出文件夹名字
  --title 图片标题
  --genename 基因名文件
	--help    usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$degfile))   { print_usage(spec)}
if ( is.null(opt$term2genefile))  { print_usage(spec) }
if ( is.null(opt$term2namefile))  { print_usage(spec) }
if ( is.null(opt$outname))  { print_usage(spec) }
if(is.null(opt$outdir)){print_usage(spec)}
if(is.null(opt$title)){print_usage(spec)}
if ( is.null(opt$genename))  { print_usage(spec) }

library(org.Ddemo.eg.db,lib.loc=opt$db)
library(clusterProfiler)
library(tidyverse)

DEG <- read.table(opt$degfile,sep = '\t',header = F, quote = '')[,1]
term2gene <- read.table(opt$term2genefile,sep = '\t',header = F, quote = '')
term2name <- read.table(opt$term2namefile,sep = '\t',header = F, quote = '')
Species_db <- "org.Ddemo.eg.db"
outname <- opt$outname
genename <- read.table(opt$genename,sep = '\t',header = T, quote = '')
if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}
setwd(opt$outdir)
###GO富集
dir.create('GO_result', showWarnings = FALSE)
GO_enrich <- enrichGO(    gene         = DEG,
                          OrgDb         = Species_db,
                          keyType       = 'GID',
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 1,
                          qvalueCutoff  = 1,
			  minGSSize = 0,
			  maxGSSize = Inf,
                          readable      = FALSE)
GO_enrich_0.05 <- enrichGO(gene         = DEG,
                           OrgDb         = Species_db,
                           keyType       = 'GID',
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
			   minGSSize = 0,
			   maxGSSize = Inf,
                           readable      = FALSE)

GO_enrich_df <- as.data.frame(GO_enrich)

#做富集气泡图
if(nrow(GO_enrich_df) > 20){GO_enrich_df <- GO_enrich_df[1:20,]}
GO_enrich_df <- GO_enrich_df[order(GO_enrich_df$ONTOLOGY),]
GO_enrich_df$Description <- factor(GO_enrich_df$Description, levels = GO_enrich_df$Description )
GO_enrich_df <- GO_enrich_df %>%
  mutate(GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio)))
p1 <- ggplot(GO_enrich_df, aes(x=GeneRatio_result,y=reorder(stringr::str_wrap(Description, 60), GeneRatio_result))) +
  geom_point(aes(size=Count, color=p.adjust, shape=ONTOLOGY)) +
  coord_cartesian(clip="off") +
  scale_color_gradient(low="red", high="yellow",name="p.adjust") +
  theme_bw() +
  scale_size_continuous(name="Count",range=c(2,10))+
  labs(x="GeneRatio", y="Terms")+
  theme(
    plot.title = element_text(hjust = 0.5,size = 18),
    axis.text.y = element_text(size = 18),  
    axis.text.x = element_text(size = 18),  
    strip.text.x = element_text(size = 18), 
    strip.text.y = element_text(size = 18), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 18) 
  )
p1=p1+ggtitle(paste("GO enrichment for",opt$title,sep=" ")) +theme(plot.title = element_text(hjust = 0.5))
ggsave(p1, file = file.path(paste0('GO_result/',outname, "_GOenrichment.png")), dpi = 600, height = 16, width = 15)
ggsave(p1, file = file.path(paste0('GO_result/',outname, "_GOenrichment.pdf")), dpi = 600, height = 16, width = 15)


#### 保存表
GO_enrich_df <- as.data.frame(GO_enrich)
# 包裹 GeneRatio 和 BgRatio 列
GO_enrich_df$GeneRatio <- paste0("(", GO_enrich_df$GeneRatio, ")")
GO_enrich_df$BgRatio <- paste0("(", GO_enrich_df$BgRatio, ")")
# 添加 gene_name 列，根据多个 Transcript_id 匹配 gene_name
GO_enrich_df$gene_name <- sapply(GO_enrich_df$geneID, function(transcript_ids) {
    # 拆分多个 Transcript_id
    transcript_ids_list <- unlist(strsplit(transcript_ids, "/"))
    # 查找每个 Transcript_id 对应的 GeneName
    gene_names <- genename$GeneName[match(transcript_ids_list, genename$transcript_id)]
    # 如果有未找到的 Transcript_id，则返回 NA
    gene_names[is.na(gene_names)] <- NA
    # 将找到的 GeneName 合并成一个字符串，以逗号分隔
    gene_names_unique <- unique(gene_names)
    paste(gene_names_unique, collapse = "/")
})
# 更改列名 geneID 为 Transcript_id
colnames(GO_enrich_df)[which(colnames(GO_enrich_df) == "geneID")] <- "Transcript_id"
# 按照 p.adjust 排序并选择需要的列
GO_enrich_df <- GO_enrich_df[order(GO_enrich_df$p.adjust), ] %>%
    select("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Transcript_id", "Count","gene_name")
write.table(GO_enrich_df,paste0('GO_result/',outname, "_GOenrichment.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)

GO_enrich_0.05_df <- as.data.frame(GO_enrich_0.05)
# 包裹 GeneRatio 和 BgRatio 列
if (nrow(GO_enrich_0.05_df) > 1) {
  GO_enrich_0.05_df$GeneRatio <- paste0("(", GO_enrich_0.05_df$GeneRatio, ")")
  GO_enrich_0.05_df$BgRatio <- paste0("(", GO_enrich_0.05_df$BgRatio, ")")
  # 添加 gene_name 列，根据多个 Transcript_id 匹配 gene_name
  GO_enrich_0.05_df$gene_name <- sapply(GO_enrich_0.05_df$geneID, function(transcript_ids) {
    # 拆分多个 Transcript_id
    transcript_ids_list <- unlist(strsplit(transcript_ids, "/"))
    # 查找每个 Transcript_id 对应的 GeneName
    gene_names <- genename$GeneName[match(transcript_ids_list, genename$transcript_id)]
    # 如果有未找到的 Transcript_id，则返回 NA
    gene_names[is.na(gene_names)] <- NA
    # 将找到的 GeneName 合并成一个字符串，以逗号分隔
    gene_names_unique <- unique(gene_names)
    paste(gene_names_unique, collapse = "/")
  })
} else {
    GO_enrich_0.05_df$gene_name <- character(0)  # 添加空列
    print("GO_enrich_0.05_df 无数据，跳过包裹操作")
}
# 更改列名 geneID 为 Transcript_id
colnames(GO_enrich_0.05_df)[which(colnames(GO_enrich_0.05_df) == "geneID")] <- "Transcript_id"
# 按照 p.adjust 排序并选择需要的列
GO_enrich_0.05_df <- GO_enrich_0.05_df[order(GO_enrich_0.05_df$p.adjust),] %>%
    select("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Transcript_id", "Count","gene_name") #按照p.adjust排序
write.table(GO_enrich_0.05_df,paste0('GO_result/',outname, "_GOenrichment_0.05.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)





###KEGG富集
if(sum(DEG %in% term2gene$V2)==0){
  print("无KEGG前景值，跳过KEGG富集")
}else{
  dir.create('KEGG_result', showWarnings = FALSE)
  KEGG_enrich <- enricher(gene = DEG,
                        TERM2GENE = term2gene,
                        TERM2NAME = term2name,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
			minGSSize = 0,
			maxGSSize = Inf)
  KEGG_enrich_0.05 <- enricher(gene = DEG,
                        TERM2GENE = term2gene,
                        TERM2NAME = term2name,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 1,
			minGSSize = 0,
			maxGSSize = Inf)
  KEGG_enrich_df <- as.data.frame(KEGG_enrich)
  #做富集气泡图
  if(nrow(KEGG_enrich_df) > 20){KEGG_enrich_df <- KEGG_enrich_df[1:20,]}
  KEGG_enrich_df$Description <- factor(KEGG_enrich_df$Description, levels = KEGG_enrich_df$Description )
  KEGG_enrich_df <- KEGG_enrich_df %>%
  mutate(GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio)))
  p2 <- ggplot(KEGG_enrich_df, aes(x=GeneRatio_result,y=reorder(stringr::str_wrap(Description, 60), GeneRatio_result))) +
  geom_point(aes(size=Count, color=p.adjust)) +
  coord_cartesian(clip="off") +
  scale_color_gradient(low="red", high="yellow",name="p.adjust") +
  theme_bw() +
  scale_size_continuous(name="Count",range=c(2,10))+
  labs(x="GeneRatio", y="Terms")+
  theme(
    plot.title = element_text(hjust = 0.5,size = 18),
    axis.text.y = element_text(size = 18),  
    axis.text.x = element_text(size = 18),  
    strip.text.x = element_text(size = 18), 
    strip.text.y = element_text(size = 18), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 18) 
  )
  p2=p2+ggtitle(paste("KEGG enrichment for",opt$title,sep=" ")) +theme(plot.title = element_text(hjust = 0.5))
  ggsave(p2, file = file.path(paste0('KEGG_result/',outname, "_KEGGenrichment.png")), dpi = 600, height = 16, width = 15)
  ggsave(p2, file = file.path(paste0('KEGG_result/',outname, "_KEGGenrichment.pdf")), dpi = 600, height = 16, width = 15)
}

  ## 保存表 
  KEGG_enrich_df <- as.data.frame(KEGG_enrich)
  # 包裹 GeneRatio 和 BgRatio 列
  KEGG_enrich_df$GeneRatio <- paste0("(", KEGG_enrich_df$GeneRatio, ")")
  KEGG_enrich_df$BgRatio <- paste0("(", KEGG_enrich_df$BgRatio, ")")
  # 添加 gene_name 列，根据多个 Transcript_id 匹配 gene_name
  KEGG_enrich_df$gene_name <- sapply(KEGG_enrich_df$geneID, function(transcript_ids) {
    # 拆分多个 Transcript_id
    transcript_ids_list <- unlist(strsplit(transcript_ids, "/"))
    # 查找每个 Transcript_id 对应的 GeneName
    gene_names <- genename$GeneName[match(transcript_ids_list, genename$transcript_id)]
    # 如果有未找到的 Transcript_id，则返回 NA
    gene_names[is.na(gene_names)] <- NA
    # 将找到的 GeneName 合并成一个字符串，以逗号分隔
    gene_names_unique <- unique(gene_names)
    paste(gene_names_unique, collapse = "/")
  })
  # 更改列名 geneID 为 Transcript_id
  colnames(KEGG_enrich_df)[which(colnames(KEGG_enrich_df) == "geneID")] <- "Transcript_id"
  KEGG_enrich_df <- KEGG_enrich_df[order(KEGG_enrich_df$p.adjust),] %>%
    select("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Transcript_id", "Count","gene_name") #按照p.adjust排序
  write.table(KEGG_enrich_df,paste0('KEGG_result/',outname, "_KEGGenrichment.xls"),sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)



  KEGG_enrich_0.05_df <- as.data.frame(KEGG_enrich_0.05)
  # 检查 KEGG_enrich_0.05_df 是否包含数据
  if (nrow(KEGG_enrich_0.05_df) > 0) {
      # 包裹 GeneRatio 和 BgRatio 列
      KEGG_enrich_0.05_df$GeneRatio <- paste0("(", KEGG_enrich_0.05_df$GeneRatio, ")")
      KEGG_enrich_0.05_df$BgRatio <- paste0("(", KEGG_enrich_0.05_df$BgRatio, ")")
      # 添加 gene_name 列，根据多个 Transcript_id 匹配 gene_name
      KEGG_enrich_0.05_df$gene_name <- sapply(KEGG_enrich_0.05_df$geneID, function(transcript_ids) {
          # 拆分多个 Transcript_id
          transcript_ids_list <- unlist(strsplit(transcript_ids, "/"))
          # 查找每个 Transcript_id 对应的 GeneName
          gene_names <- genename$GeneName[match(transcript_ids_list, genename$transcript_id)]
          # 如果有未找到的 Transcript_id，则返回 NA
          gene_names[is.na(gene_names)] <- NA
          # 将找到的 GeneName 合并成一个字符串，以逗号分隔
          gene_names_unique <- unique(gene_names)
          paste(gene_names_unique, collapse = "/")
      })
  } else {
      print("KEGG_enrich_0.05_df 无数据，跳过包裹操作")
      # KEGG_enrich_0.05_df$gene_name <- NA  # 空数据框时添加空列
      KEGG_enrich_0.05_df$gene_name <- character(0)  # 添加空列
  }
  # 更改列名 geneID 为 Transcript_id
  if ("geneID" %in% colnames(KEGG_enrich_0.05_df)) {
      colnames(KEGG_enrich_0.05_df)[which(colnames(KEGG_enrich_0.05_df) == "geneID")] <- "Transcript_id"
  } else {
      print("geneID 列不存在，无法更改列名")
  }
  # 按照 p.adjust 排序，并选择需要的列
  KEGG_enrich_0.05_df <- KEGG_enrich_0.05_df[order(KEGG_enrich_0.05_df$p.adjust),] %>%
          select("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Transcript_id", "Count", "gene_name")
  # 保存结果
  write.table(KEGG_enrich_0.05_df, paste0('KEGG_result/', outname, "_KEGGenrichment_0.05.xls"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)






escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()



