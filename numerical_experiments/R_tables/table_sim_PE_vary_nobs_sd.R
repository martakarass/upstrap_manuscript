
#' Compute table summary of PE 

rm(list = ls())
library(tidyverse)
library(here)
library(stargazer)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# READ AND PROCESS DATA 

func_get_table <- function(dat){
  name_tmp <- c("upstrap_power", "pe_ups_true")
  tbl <- 
    dat %>%
    filter(name %in% name_tmp) %>%
    mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
    mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
    mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
    select(N_obs, N_tar, sd_sigma, name, value_f) %>%
    pivot_wider(names_from = name, values_from = value_f) %>%
    arrange(N_obs, N_tar) %>%
    mutate(sd_sigma = as.character(sd_sigma)) %>%
    # fix 2021-12-11 to have "sample size" denote all units (not: units per treatment arm)
    mutate(N_obs = 2 * N_obs) %>% 
    mutate(N_tar = 2 * N_tar) %>% 
    as.data.frame() %>%
    select(N_obs, N_tar, sd_sigma, upstrap_power, pe_ups_true)
  tbl[is.na(tbl)] <- "-"
  tbl
}


# ------------------------------------------------------------------------------
# PROBLEM 10

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-twosample_ttest_vary_nobs_sd_agg.rds")
dat <- readRDS(fpath_tmp)
sort(unique(dat$sd_sigma))
dat <- dat %>% filter(sd_sigma == sort(unique(dat$sd_sigma))[2])
dat_10_tbl <- func_get_table(dat)
dat_10_tbl

  
# ------------------------------------------------------------------------------
# PROBLEM 11

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lm_testcoef_vary_nobs_sd_agg.rds")
sort(unique(dat$sd_sigma))
dat <- readRDS(fpath_tmp)
dat <- dat %>% filter(sd_sigma == sort(unique(dat$sd_sigma))[2])
dat_11_tbl <- func_get_table(dat)
dat_11_tbl


# ------------------------------------------------------------------------------
# PROBLEM 12

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lmm_testcoef_vary_nobs_sd_agg.rds")
sort(unique(dat$sd_sigma))
dat <- readRDS(fpath_tmp)
dat <- dat %>% filter(sd_sigma == sort(unique(dat$sd_sigma))[2])
dat_12_tbl <- func_get_table(dat)
dat_12_tbl


# ------------------------------------------------------------------------------

library(stargazer)

stargazer(dat_10_tbl, summary = FALSE)
stargazer(dat_11_tbl, summary = FALSE)
stargazer(dat_12_tbl, summary = FALSE)




