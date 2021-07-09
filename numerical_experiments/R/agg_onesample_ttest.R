
#' This script aggregate power of rejecting H0 in one-sample t-test problem. 
#' 
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave agg_onesample_ttest.R -N JOB_agg_onesample

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-onesample_ttest_raw")

# read and combine data 
fnames <- list.files(res_fdir_raw, full.names = TRUE)
fnames <- fnames[grepl("arrayjob", fnames)]
dat_l <- lapply(fnames, readRDS)
dat <- do.call("rbind", dat_l)
dat <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))

length(fnames); length(unique(dat$arrayjob_idx))
dim(dat)
head(dat)
table(dat$name)
table(dat$eff_tar, useNA = "always")

# aggregate data
dat_agg <- 
  dat %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value)
  ) %>%
  ungroup()

dat_diff <- 
  dat %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(upstrap_powerttest_power_diff = upstrap_power - powerttest_power) %>%
  select(-c(upstrap_power, powerttest_power)) %>%
  pivot_longer(cols = upstrap_powerttest_power_diff)
dat_agg_diff <- 
  dat_diff %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value)
  ) %>%
  ungroup()

dat_out <- rbind(dat_agg, dat_agg_diff)
dim(dat_out)
table(dat_out$name)

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-onesample_ttest_agg.rds")
saveRDS(dat_out, fpath_tmp)


