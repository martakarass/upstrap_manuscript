#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave adhoc_test_glmm_testcoef.R -t 1-50 -tc 50 -N JOB_adhoc_glmm
#' ls -l-d *adhoc*


dir_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-04-glmm_adhoc_test")

fnames <- list.files(dir_out, full.names = TRUE)
dat <- do.call("rbind", lapply(fnames, readRDS)) 

dat
