
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
  N_obs_tmp <- 50
  name_tmp <- c("upstrap_power", "comparator_power", "pe_ups_true", "pe_cmp_true", "pe_ups_cmp")
  tbl <- 
    dat %>%
    filter(N_obs == N_obs_tmp, name %in% name_tmp) %>%
    mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
    mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
    mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
    select(eff_tar, N_tar, name, value_f) %>%
    pivot_wider(names_from = name, values_from = value_f) %>%
    arrange(eff_tar, N_tar) %>%
    # fix 2021-12-11 to have "sample size" denote all units (not: units per treatment arm)
    mutate(N_tar = 2 * N_tar) %>% 
    as.data.frame() %>%
    select(eff_tar, N_tar, upstrap_power, comparator_power, pe_ups_cmp, pe_ups_true, pe_cmp_true)
  tbl[is.na(tbl)] <- "-"
  tbl
}


# ------------------------------------------------------------------------------
# PROBLEM 1

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-onesample_ttest_agg.rds")
dat <- readRDS(fpath_tmp)
dat_1_tbl <- func_get_table(dat)

dat_1_tbl

  
# ------------------------------------------------------------------------------
# PROBLEM 2

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-twosample_ttest_agg.rds")
dat <- readRDS(fpath_tmp)
table(dat$name)
dat_2_tbl <- func_get_table(dat)

dat_2_tbl


# ------------------------------------------------------------------------------
# PROBLEM 3

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lm_testcoef_agg.rds")
dat <- readRDS(fpath_tmp)
table(dat$name)
dat_3_tbl <- func_get_table(dat)

dat_3_tbl


# ------------------------------------------------------------------------------
# PROBLEM 4

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glm_testcoef_agg.rds")
dat <- readRDS(fpath_tmp)
table(dat$name)
dat_4_tbl <- func_get_table(dat)

dat_4_tbl


# ------------------------------------------------------------------------------
# PROBLEM 5

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lmm_testcoef_agg.rds")
dat <- readRDS(fpath_tmp)
table(dat$name)
dat_5_tbl <- func_get_table(dat)

dat_5_tbl



# ------------------------------------------------------------------------------
# PROBLEM 6

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glmm_testcoef_agg.rds")
dat <- readRDS(fpath_tmp)
table(dat$name)
dat_6_tbl <- func_get_table(dat)

dat_6_tbl


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(stargazer)

stargazer(dat_1_tbl, summary = FALSE)
stargazer(dat_2_tbl, summary = FALSE)
stargazer(dat_3_tbl, summary = FALSE)
stargazer(dat_4_tbl, summary = FALSE)
stargazer(dat_5_tbl, summary = FALSE)
stargazer(dat_6_tbl, summary = FALSE)




