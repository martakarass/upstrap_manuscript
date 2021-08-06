#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave adhoc_test_glmm_testcoef.R -t 1-50 -tc 50 -N JOB_adhoc_glmm
#' ls -l-d *adhoc*

rm(list = ls())
library(tidyverse)
dir_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-04-glmm_adhoc_test")

fnames <- list.files(dir_out, full.names = TRUE)
length(fnames)

dd <- readRDS(fnames[1])
head(dd)
"lanuch" %in% colnames(dd)

dat <- do.call("rbind", lapply(fnames, function(fname){
  dd <- readRDS(fname)
  if ("lanuch" %in% colnames(dd)){
    return(dd)
  } else {
    return(NULL)
  }
})) 
dat <- dat %>% filter(lanuch == 3, !is.na(power_upstrap_v1)) 

head(dat)
dim(dat)


dim(dat); dat %>% summarise_all(mean, rm.na = TRUE)

# ---> launch 0
# 483   5
#   result_glmm power_upstrap power_upstrap_nore power_simr arrayjob_idx
# 1   0.4720497     0.6045507          0.6943892   0.542058     24.65839

dim(dat); dat %>% summarise_all(mean, rm.na = TRUE)

# ---> launch 1
# [1] 415   7
# result_glmm power_upstrap_v1 power_upstrap_v2 power_upstrap_v3 power_simr
# 1   0.4819277        0.6049157               NA         0.542612  0.5401639
# arrayjob_idx lanuch
# 1     24.02169      2

dim(dat); dat %>% summarise_all(mean, rm.na = TRUE)



