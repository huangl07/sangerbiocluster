#!/usr/bin/env Rscript

suppressMessages({
library(dplyr)
library(rmarkdown)
library(kableExtra)
library(rmdformats)
library(bookdown)
library(latex2exp)
library(argparser)})
argv <- arg_parser('插入位点检测报告')
argv <- add_argument(argv,"--rmd", help="the rmd file")
argv <- add_argument(argv,"--format", help="the output format")
argv <- add_argument(argv,"--dev", flag=TRUE, help="dev mode")
argv <- add_argument(argv,"--outfile", help="the output file")

argv <- parse_args(argv)
format <- argv$format
outfile <- argv$outfile
rmdfile <- argv$rmd

options(bitmapType='cairo-png')

if (format=='html') {
    render(input=rmdfile,
        output_format='bookdown::html_document2',
        output_file=outfile,
        params = list(
            docfmt = "html",
            devmode = argv$dev
    )
)} else if (format=='pdf') {
    render(input=rmdfile,
        output_format='bookdown::pdf_document2',
        output_file=outfile,
        params = list(
            docfmt = "pdf",
            devmode = argv$dev
    )
)} else {
    print("format error!")
    q(1)
}

