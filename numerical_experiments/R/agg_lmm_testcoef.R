
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave agg_lm_testcoef.R -N JOB_agg_lm_testcoef

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-lm_testcoef_raw")

fnames <- list.files(res_fdir_raw, full.names = TRUE)
length(fnames)

dat <- lapply(fnames, readRDS)
dat <- do.call("rbind", dat)
dim(dat)
head(dat)
table(dat$name)
table(dat$eff_tar, useNA = "always")

dat_agg <- 
  dat %>%
  # mutate(eff_tar = ifelse(is.na(eff_tar)))
  group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value)
  )

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-lm_testcoef_agg.rds")
saveRDS(dat_agg, fpath_tmp)
