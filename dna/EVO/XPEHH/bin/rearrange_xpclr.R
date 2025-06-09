library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-x", "--xpclr"), type = "character", action = "store", help = "Input xpclr file"),
  make_option(c("-t", "--thread"), type = "integer", action = "store", default = 8, help = "CPUS [default %default]")
)
opt = parse_args(OptionParser(option_list = option_list,
                              usage = "Usage: %prog [options] \nDescription: This Script is used to rearrange xpclr result!"))

times<-Sys.time()

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
suppressMessages(library(furrr))
suppressMessages(library(progressr))
suppressMessages(library(magrittr))

write_result_fun <-function(xpclr_abstract){
  xpclr_abstract <- xpclr_abstract %>% 
  write.table(file=paste(xpclr_abstract$pop[1], "xpclr_abstract",sep="."),row.names=FALSE, quote = FALSE, sep = "\t", xpclr_abstract)
}


###批量读取tajima文件, pi文件, fst文件, dxy文件
files <- fs::dir_ls(opt$xpclr, recurse = TRUE, glob = "*-*.xpclr")
xpclr <- map_dfr(files, fread, .id = "pop", sep="\t", header=TRUE, integer64 = "numeric")

xpclr$pop <- xpclr$pop %>% 
  gsub(".*/","", .) %>% 
  str_replace("^(.+)-LG.*", "\\1")
xpclr_abstract <- xpclr %>%  
  select(pop, id, chrom, pos_start, pos_stop,xpclr, xpclr_norm)
xpclr_abstract_list <- xpclr_abstract %>% 
  group_split(pop)

plan(multisession, workers = opt$thread) #计划使用8个CPU
with_progress({
  pb <- progressor(steps = length(xpclr_abstract_list)) #进度条，R>4.0.0
  xpclr_abstract_list %>% 
    future_walk(write_result_fun)
})

escape_time <- Sys.time()-times
print("Done!")
print(escape_time)
