#!/usr/bin/env Rscript

suppressMessages({
library(dplyr)
library(rmarkdown)
library(kableExtra)
library(rmdformats)
library(bookdown)
library(latex2exp)
library(argparser)})
argv <- arg_parser("BSA生成报告")
argv <- add_argument(argv, "--rmd", help="the rmd file")
argv <- add_argument(argv, "--format", help="the output format")
argv <- add_argument(argv, "--outfile", help="the output file")
argv <- add_argument(argv, "--population", help="群体")
argv <- add_argument(argv, "--primer_design", help = "是否有引物设计")
argv <- add_argument(argv, "--deepbsa", help = "是否有deepbsa")
argv <- add_argument(argv, "--mutmap", help = "是否为mutmap")
argv <- parse_args(argv)
format <- argv$format
outfile <- argv$outfile
rmdfile <- argv$rmd
population <- argv$population
primer_design <- argv$primer_design

options(bitmapType = 'cairo-png')

if (format == "html") {
    render(input = rmdfile,
        output_format = "bookdown::html_document2",
        output_file = outfile,
        params = list(
            docfmt = "html",
            population = population,
            primer_design = primer_design,
            mutmap = argv$mutmap,
            deepbsa = argv$deepbsa
    )
)} else if (format=='pdf') {
    render(input=rmdfile,
        output_format='bookdown::pdf_document2',
        output_file=outfile,
        params = list(
            docfmt = "pdf",
            population = population,
            primer_design = primer_design,
            mutmap = argv$mutmap,
            deepbsa = argv$deepbsa
    )
)} else {
    print("format error!")
    q(1)
}
