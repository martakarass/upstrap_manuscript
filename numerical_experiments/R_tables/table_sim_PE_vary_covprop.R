
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
    select(cov_prop, N_tar, name, value_f) %>%
    pivot_wider(names_from = name, values_from = value_f) %>%
    arrange(cov_prop, N_tar) %>%
    mutate(cov_prop = as.character(cov_prop)) %>%
    # fix 2021-12-11 to have "sample size" denote all units (not: units per treatment arm)
    mutate(N_tar = 2 * N_tar) %>% 
    as.data.frame() %>%
    select(cov_prop, N_tar, upstrap_power, pe_ups_true)
  tbl[is.na(tbl)] <- "-"
  tbl
}


# ------------------------------------------------------------------------------
# PROBLEM 7

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-twosample_ttest_vary_covprop_agg.rds")
dat <- readRDS(fpath_tmp)
dat_7_tbl <- func_get_table(dat)
dat_7_tbl

  
# ------------------------------------------------------------------------------
# PROBLEM 8

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lm_testcoef_vary_covprop_agg.rds")
dat <- readRDS(fpath_tmp)
dat_8_tbl <- func_get_table(dat)
dat_8_tbl


# ------------------------------------------------------------------------------
# PROBLEM 9

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lmm_testcoef_vary_covprop_agg.rds")
dat <- readRDS(fpath_tmp)
dat_9_tbl <- func_get_table(dat)
dat_9_tbl


# ------------------------------------------------------------------------------

library(stargazer)

stargazer(dat_7_tbl, summary = FALSE)
stargazer(dat_8_tbl, summary = FALSE)
stargazer(dat_9_tbl, summary = FALSE)




