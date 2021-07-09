
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

dat_l <- lapply(fnames, readRDS)
dat <- do.call("rbind", dat_l)
dat <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))

table(dat$name)
table(dat$name, dat$eff_tar)

# aggregate power 
dat_agg <- 
  dat %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value)
  ) %>%
  ungroup()

table(dat_agg$name, dat_agg$eff_tar)


# aggregate power difference
dat_diff <- 
  dat %>%
  filter(name %in% c("upstrap_power", "simr_power")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(upstrap_simr_power_diff = upstrap_power - simr_power) %>%
  select(-c(upstrap_power, simr_power)) %>%
  pivot_longer(cols = upstrap_simr_power_diff)
dat_agg_diff <- 
  dat_diff %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value)
  ) %>%
  ungroup()

table(dat_agg_diff$name, dat_agg_diff$eff_tar)


dat_out <- rbind(dat_agg, dat_agg_diff)
dim(dat_out)
table(dat_out$name)

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-lm_testcoef_agg.rds")
saveRDS(dat_out, fpath_tmp)
