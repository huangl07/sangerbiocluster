library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-p", "--pi"), type = "character", action = "store",default=NULL, help = "Input pi file"),
  make_option(c("-f", "--fst"), type = "character", action = "store",default=NULL, help = "Input fst file"),
  make_option(c("-d", "--dxy"), type = "character", action = "store",default=NULL, help = "Input dxy file"),
  make_option(c("-t", "--tajimaD"), type = "character", action = "store",default=NULL, help = "Input tajimaD dir"),
  make_option(c("-m", "--mult"), type = "integer", action = "store", default = 8, help = "CPUS [default %default]")
)
times<-Sys.time()

opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"))
if(is.null(opt$pi)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$fst)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$dxy)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$tajimaD)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}



suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
suppressMessages(library(furrr))
suppressMessages(library(magrittr))

merge_fst_pi_tajimaD_fun <-function(fst, pi_abstract, tajimaD_abstract){
  pi1 <- pi_abstract %>% filter(pop == fst$pop1[1]) #pop1
  pi2 <- pi_abstract %>% filter(pop == fst$pop2[1]) #pop2
  tajimaD1 <- tajimaD_abstract %>% filter(pop == fst$pop1[1]) #tajimaD1
  tajimaD2 <- tajimaD_abstract %>% filter(pop == fst$pop2[1]) #tajimaD2
  ##合并fst,pi1,pi2
  merge_result <- fst %>% 
    merge(pi1, 
          all.x = TRUE,
          sort = FALSE, 
          by.x = c("pop1", "chromosome", "window_pos_1", "window_pos_2"), 
          by.y = c("pop", "chromosome", "window_pos_1", "window_pos_2")) %>% 
    rename(pi1 = avg_pi) %>% 
    merge(pi2, 
          all.x = TRUE,
          sort = FALSE, 
          by.x = c("pop2", "chromosome", "window_pos_1", "window_pos_2"), 
          by.y = c("pop", "chromosome", "window_pos_1", "window_pos_2")) %>% 
    rename(pi2 = avg_pi) %>% 
    merge(tajimaD1, 
          all.x = TRUE,
          sort = FALSE, 
          by.x = c("pop1", "chromosome", "window_pos_1", "window_pos_2"), 
          by.y = c("pop", "chromosome", "window_pos_1", "window_pos_2")) %>% 
    rename(tajimaD1 = tajimaD) %>% 
    merge(tajimaD2, 
          all.x = TRUE,
          sort = FALSE, 
          by.x = c("pop2", "chromosome", "window_pos_1", "window_pos_2"), 
          by.y = c("pop", "chromosome", "window_pos_1", "window_pos_2")) %>% 
    rename(tajimaD2 = tajimaD) %>% 
    select(pop1, pop2, chromosome, window_pos_1, window_pos_2, pi1, pi2, avg_wc_fst, tajimaD1, tajimaD2, avg_dxy)
    merge_result$avg_wc_fst[merge_result$avg_wc_fst < 0] = 0
    merge_result %>% 
      mutate(z_avg_wc_fst = scale(avg_wc_fst)) %T>% 
    write.table(file=paste(fst$pop1[1], fst$pop2[1], "detail",sep="."), row.names=FALSE, quote = FALSE, sep = "\t") %>% 
    return()
}

###批量读取tajima文件, pi文件, fst文件, dxy文件
files <- fs::dir_ls(opt$tajimaD, recurse = TRUE, glob = "*.Tajima.D")
tajimaD <- map_dfr(files, fread, .id = "pop", sep="\t", header=TRUE, integer64 = "numeric")
pi <- fread(opt$pi ,sep="\t",header=TRUE, integer64 = "numeric")
fst <- fread(opt$fst ,sep="\t",header=TRUE, integer64 = "numeric")
dxy <- fread(opt$dxy, sep="\t", header=TRUE, integer64 = "numeric")

###整理输入文件，fs文件为基准文件
pi_abstract <- pi %>% 
  select(pop, chromosome, window_pos_1, window_pos_2, avg_pi) %>% 
  filter(avg_pi != "NA", avg_pi != "NaN")
dxy_abstract <- dxy %>% 
  select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_dxy) %>% 
  filter(avg_dxy != "NA", avg_dxy != "NaN")
##正则提取pop
tajimaD$pop <- tajimaD$pop %>% 
  gsub(".*/","", .) %>% 
  str_replace("^(.+)\\.Tajima\\.D", "\\1")
tajimaD_abstract <- tajimaD %>%
  mutate(window_pos_2 = BIN_START + 10000) %>% 
  mutate(window_pos_1 = BIN_START + 1) %>% 
  rename(chromosome = CHROM, tajimaD = TajimaD) %>% 
  select(pop, chromosome, window_pos_2, window_pos_1, tajimaD) %>% 
  filter(tajimaD != "NA", tajimaD != "NaN")
fst_abstract <- fst %>%  
  select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_wc_fst) %>% 
  filter(avg_wc_fst != "NA", avg_wc_fst != "NaN")

###合并文件
##合并fst和dxy
merge_result_list <- fst_abstract %>%
  merge(dxy_abstract, 
        all.x = TRUE, 
        sort = FALSE,
        by =c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")) %>% 
  group_split(pop1, pop2)

plan(multisession, workers = opt$mult) #计划使用8个CPU
merge_result_list %>% 
  future_walk(merge_fst_pi_tajimaD_fun, pi_abstract, tajimaD_abstract)


escape_time <- Sys.time()-times
print("Done!")
print(escape_time)
