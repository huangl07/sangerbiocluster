#!/usr/bin/env Rscript
library(argparser)
library(ggplot2)
times <- Sys.time()
argv <- arg_parser("plot_phenolyzer.R")
argv <- add_argument(argv, "--infile", help = "phenolyzer result file")
argv <- add_argument(argv, "--outfile", help = "the output file")
argv <- add_argument(argv, "--dfile", help="disease file")
argv <- parse_args(argv)

disease <- readLines(argv$dfile)[1]
df <- read.delim(argv$infile, sep = "\t", header = TRUE)
if (nrow(df) > 0) {
    names(df)[1] <- "Disease"
    df$Disease <- disease
    write.table(df, paste0(argv$outfile, ".xls"), sep = "\t", row.names = F, col.names = T, quote = F)
    df$Gene <- factor(df$Gene, levels = rev(df$Gene))
    plot <- ggplot(df[seq_len(min(20, nrow(df))), ]) +
        geom_col(aes(fill = Status, y = Gene, x = Score), color = NA, orientation = "y") +
        xlim(0, 1) + labs(title=paste("disease contributory genes of", disease))+
        scale_fill_manual(name = "Gene Status", breaks = c("Predicted", "SeedGene"), values = c("orange2", "chartreuse3")) +
        theme(
            panel.border = element_rect(fill = NA, color = "black", linetype = 1),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(0.25, "cm"),
            axis.line = element_blank(),
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 15),
            axis.title.x = element_text(margin = margin(0.5, 0, 0, 0, "cm")),
            axis.title.y = element_text(margin = margin(0, 0.5, 0, 0, "cm")),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.key = element_blank(),
            legend.justification = c(1, 1),
            legend.background = element_blank()
        )
    ggsave(paste0(argv$outfile, ".pdf"), plot, dpi = 600, width = 10, height = 10)
    ggsave(paste0(argv$outfile, ".png"), plot, dpi = 600, width = 10, height = 10)
}
escaptime <- Sys.time() - times
print("Done!")
print(escaptime)
q()
