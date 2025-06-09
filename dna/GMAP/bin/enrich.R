#!/usr/bin/env Rscript
# @Last-edit Time 2022/12/1
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com
options(bitmapType = "cairo")
library(getopt)
spec <- matrix(c(
    "degfile", "d", 1, "character",
    "gson", "g", 2, "character",
    "degcol","c",2,"character",
    "outname", "o", 1, "character",
        "wgs", "w", 2, "character",
        "prefix","p",1,"character",
    "rlib", "b", 2, "character",
        "godb", "s", 2, "character",
    "type", "t", 2, "character",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec = spec) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("
Usage:
        --degfile,-d   输入基因列表, '-'表示从标准输入中读取
        --degcol,-c    输入deg对应列名
        --gson,-g   输入kegg.gson文件
        --outname,-o   输出文件夹
        --prefix,-p     输出前缀
        --godb,-s       go注释db包名, default org.Ddemo.eg.db
        --rlib,-b   输入go注释db包所在lib路径
        --type,-t   kegg or go, default both
        --wgs,-w        wgsv4工作目录，用来抓取gson和rlib，可不填
        --help,-h   usage
\n")
    q(status = 1)
}
if (!is.null(opt$help)) {
    print_usage(spec)
}
if (is.null(opt$degfile)) {
    print_usage(spec)
} else if(opt$degfile=="-"){
        opt$degfile<-"stdin"
}
if (is.null(opt$outname)) {
    opt$outname<-"."
}
if (is.null(opt$type)) {
    runkegg <- TRUE
        rungo <- TRUE
} else if (opt$type == "kegg"){
        runkegg <- TRUE
        rungo <- FALSE
} else if (opt$type == "go"){
        runkegg <- FALSE
    rungo <- TRUE
} else {
        print_usage(spec)
}
if (!is.null(opt$wgs)){
        opt$gson <- file.path(opt$wgs,"output/tmp/02.reference/orgDB/kegg.gson")
        opt$rlib <- file.path(opt$wgs,"AnnoVar/AnnoAnalysis/RLib/")
}
if (is.null(opt$gson) && runkegg) {
    print_usage(spec)
}
if (is.null(opt$godb)) {
    opt$godb <- "org.Ddemo.eg.db"
}
print(opt)
library(tidyverse)
library(clusterProfiler)
library(LaF)
if(is.null(opt$degcol)){
    DEG <- readLines(opt$degfile)
}else{
    model <- detect_dm_csv(opt$degfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", nrows=10)
    laf <- laf_open(model)
    DEG<-process_blocks(laf, function(df, result){unique(c(result,df[,opt$degcol]))}, nrows = 5000)
}
DEG <- DEG[!duplicated(DEG)]
outname <- opt$outname
if (length(DEG)==0) {
        cat("no DEG in degfiles!")
        q(status=0)
}else{
        cat("Total",length(DEG),"DEGs in degfiles.",sep=" ")
}
if (rungo) {
    ### go富集
    library(opt$godb,character.only=TRUE, lib.loc = opt$rlib)
    go_enrich <- enrichGO(
        gene = DEG,
        OrgDb = opt$godb,
        keyType = "GID",
        ont = "ALL",
        pAdjustMethod = "BH",
        maxGSSize = Inf,
        minGSSize = 0,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = FALSE
    )
    go_enrich_df <- as.data.frame(go_enrich)
    go_enrich_df <- go_enrich_df[order(go_enrich_df$p.adjust), ]
    write.table(go_enrich_df,
        file.path(
            outname,
            paste0(opt$prefix,".go.enrichment.xls")
        ),
        sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
    )
    # 做go富集气泡图
    go_enrich_df <- go_enrich_df[seq_len(min(20, nrow(go_enrich_df))), ]
    go_enrich_df$logp <- -log(go_enrich_df$p.adjust, 10)
    go_enrich_df <- go_enrich_df[order(go_enrich_df$logp,decreasing=TRUE), ]
    go_enrich_df$Description <- factor(go_enrich_df$Description,
        levels = go_enrich_df$Description
    )
    go_enrich_df <- go_enrich_df %>%
        mutate(
            GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio))
        )
    p1 <- ggplot(
        go_enrich_df,
        aes(x = GeneRatio_result, y = stringr::str_wrap(Description, 40))
    ) +
        geom_point(aes(
            size = Count, color = logp,
            shape = ONTOLOGY
        )) +
        coord_cartesian(clip = "off") +
        scale_color_gradient(
            low = "yellow", high = "red",
            name = "-log10(p.adjust)"
        ) +
        theme_bw() +
        scale_size_continuous(name = "Count", range = c(2, 12)) +
        labs(x = "GeneRatio", y = "GO Terms", title = "Go Enrichment") +
        theme(plot.title = element_text(hjust=0.5))
    ggsave(p1,
        file = file.path(
            outname,
            paste0(opt$prefix,".go.enrichment.png")
        ),
        dpi = 600, height = 10, width = 15
    )
    ggsave(p1,
        file = file.path(
            outname,
            paste0(opt$prefix,".go.enrichment.pdf")
        ),
        dpi = 600, height = 10, width = 15
    )
}
if (runkegg) {
    gson <- gson::read.gson(opt$gson)
    ### kegg富集
    if (sum(DEG %in% gson@gsid2gene[, 2]) == 0) {
        print("无KEGG前景值，跳过KEGG富集")
    } else {
        kegg_enrich <- enricher(
            gene = DEG,
            TERM2GENE = gson@gsid2gene,
            TERM2NAME = gson@gsid2name,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            maxGSSize = Inf,
            minGSSize = 0
        )
        kegg_enrich_df <- as.data.frame(kegg_enrich)
        kegg_enrich_df$logp <- -log(kegg_enrich_df$p.adjust, 10)
        kegg_enrich_df <- kegg_enrich_df[order(kegg_enrich_df$logp,decreasing=TRUE), ]
        write.table(kegg_enrich_df,
            file.path(
                outname,
                paste0(opt$prefix,".kegg.enrichment.xls")
            ),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
        )
        # 做kegg富集气泡图
        kegg_enrich_df <- kegg_enrich_df[seq_len(min(20, nrow(kegg_enrich_df))), ]
        kegg_enrich_df$Description <- factor(kegg_enrich_df$Description,
            levels = kegg_enrich_df$Description
        )
        kegg_enrich_df <- kegg_enrich_df %>%
            mutate(
                GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio))
            )
        p2 <- ggplot(
            kegg_enrich_df,
            aes(x = GeneRatio_result, y = stringr::str_wrap(Description, 40))
        ) +
            geom_point(aes(size = Count, color = logp)) +
            coord_cartesian(clip = "off") +
            scale_color_gradient(
                low = "yellow", high = "red",
                name = "-log10(p.adjust)"
            ) +
            theme_bw() +
            scale_size_continuous(name = "Count", range = c(2, 12)) +
            labs(x = "GeneRatio", y = "Pathways", title = "KEGG Enrichment") +
            theme(plot.title = element_text(hjust=0.5))
        ggsave(p2,
            file = file.path(
                outname,
                paste0(opt$prefix,".kegg.enrichment.png")
            ),
            dpi = 600, height = 10, width = 15
        )
        ggsave(p2,
            file = file.path(
                outname,
                paste0(opt$prefix,".kegg.enrichment.pdf")
            ),
            dpi = 600, height = 10, width = 15
        )
    }
}

print("Done!")

