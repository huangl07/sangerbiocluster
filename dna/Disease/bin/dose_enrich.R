#!/usr/bin/env Rscript
library(DOSE)
times <- Sys.time()
library("getopt")
spec <- matrix(c(
    "gene", "g", 1, "character",
    "out", "o", 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
gene <- read.table(opt$gene, sep = "\t", header = FALSE)
if (nrow(gene) > 0) {
    x <- enrichDGN(
        gene = gene[, 1],
        pvalueCutoff = 1, # 显著可设置0.05
        pAdjustMethod = "BH",
        minGSSize = 1,
        maxGSSize = Inf,
        qvalueCutoff = 1, # 显著可设置0.05
        readable = FALSE
    )
    df <- as.data.frame(x)
    df <- df[, -which(names(df) == "qvalue")]
    df <- df[order(df$p.adjust, decreasing = FALSE), ]
    write.table(df, paste0(opt$out, ".enrichment.xls"), sep = "\t", quote = F, row.names = F)
    df <- df[seq_len(min(20, nrow(df))), ]
    df$logp <- -log(df$p.adjust, 10)
    df <- df %>%
        mutate(
            GeneRatio_result = Count / as.numeric(sub("\\d+/", "", GeneRatio))
        )
    p2 <- ggplot(
        df,
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
        labs(x = "GeneRatio", y = "Disease", title = paste("DisGeNET Enrichment of", out_title)) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(p2,
        file = paste0(opt$out, ".enrichment.png"),
        dpi = 600, height = 10, width = 15
    )
    ggsave(p2,
        file = paste0(opt$out, ".enrichment.pdf"),
        dpi = 600, height = 10, width = 15
    )
    writeLines(df[seq_len(min(5, nrow(df))), ]$Description, paste0(opt$out, ".list"))
} else {
    cat("", file = paste0(opt$out, ".list"))
}
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
