#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# rm(list = ls())

#' Aggregate results generated on the cluster via array job 
#' (each array job saved a separate file).s

rm(list = ls())
library(data.table)
library(tidyverse)
options(scipen=999)

# pull files from the 
fdir_tmp  <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-05-01-upstrap_no_agg_results/")
fname_list <- list.files(fdir_tmp, full.names = TRUE)
dat <- lapply(fname_list, readRDS)
dat <- do.call("rbind", dat)
dim(dat)
head(dat)

dat %>% 
  group_by(time_point, zone, sample_size_M) %>%
  summarise(
    power = mean(power),
    eval_time_secs = median(eval_time_secs),
    B_boot_success = mean(B_boot_success)
  ) %>%
  as.data.frame()

# save aggregated file 
fpath_tmp <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-05-01-upstrap_no_agg_results.rds")
saveRDS(dat, fpath_tmp)

# get /users/mkaras/_PROJECTS/upstrap_manuscript/use_case_examples/bausema2016/results/2021-05-01-upstrap_no_agg_results.rds /Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/use_case_examples/bausema2016/results/2021-05-01-upstrap_no_agg_results.rds
