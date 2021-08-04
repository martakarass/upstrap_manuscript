
#' Compute table summary of PE 

rm(list = ls())
library(tidyverse)
library(here)
library(stargazer)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# READ AND PROCESS DATA 

N_tar_NN <- 6

func_table_pwr <- function(dat){
  dat %>%
    filter(N_tar %in% N_tar_vals_sub) %>%
    filter(name == "upstrap_power") %>%
    mutate(value_mean_f   = sprintf("%.3f", value_mean)) %>%
    mutate(value_sd_f     = sprintf("%.3f", value_sd)) %>%
    mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
    select(N_obs, N_tar, eff_tar, value_f) %>%
    pivot_wider(names_from = eff_tar, values_from = value_f) 
}

func_table_pe <- function(dat){
  dat %>%
    filter(N_tar %in% N_tar_vals_sub) %>%
    filter(name == "ups_cmp_power_pe") %>%
    mutate(value_mean_f   = sprintf("%.3f", value_mean)) %>%
    mutate(value_sd_f     = sprintf("%.3f", value_sd)) %>%
    mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
    select(N_obs, N_tar, eff_tar, value_f) %>%
    pivot_wider(names_from = eff_tar, values_from = value_f) 
}


# ------------------------------------------------------------------------------
# PROBLEM 1

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-onesample_ttest_agg.rds")

dat <- readRDS(fpath_tmp)
table(dat$name)
N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
N_tar_vals_sub[which_repl] <- N_obs_val

dat_1_pwr <- func_table_pwr(dat)
dat_1_pe <- func_table_pe(dat)

dat_1_pwr
dat_1_pe


# ------------------------------------------------------------------------------
# PROBLEM 2

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-twosample_ttest_agg.rds")

dat <- readRDS(fpath_tmp)
table(dat$name)

N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
N_tar_vals_sub[which_repl] <- N_obs_val

dat_2_pwr <- func_table_pwr(dat)
dat_2_pe <- func_table_pe(dat)

dat_2_pwr
dat_2_pe


# ------------------------------------------------------------------------------
# PROBLEM 3

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-lm_testcoef_agg.rds")

dat <- readRDS(fpath_tmp)
head(dat)
table(dat$cnt)
table(dat$name)
table(dat$eff_tar)
sort(unique((dat$N_tar)))

N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
N_tar_vals_sub[which_repl] <- N_obs_val

dat_3_pwr <- func_table_pwr(dat)
dat_3_pe <- func_table_pe(dat)

dat_3_pwr
dat_3_pe


# ------------------------------------------------------------------------------
# PROBLEM 4

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-02-glm_testcoef_agg.rds")

dat <- readRDS(fpath_tmp)
head(dat)
table(dat$cnt)
table(dat$name)
table(dat$eff_tar)
sort(unique((dat$N_tar)))

N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
N_tar_vals_sub[which_repl] <- N_obs_val

dat_4_pwr <- func_table_pwr(dat)
dat_4_pe <- func_table_pe(dat)

dat_4_pwr
dat_4_pe



# ------------------------------------------------------------------------------
# PROBLEM 5

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-08-lmm_testcoef_agg.rds")

dat <- readRDS(fpath_tmp)
head(dat)
table(dat$cnt)
table(dat$name)
table(dat$eff_tar)
sort(unique((dat$N_tar)))

N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
repl_with <- N_obs_val
if (!(repl_with %in% N_tar_vals)){
  repl_with <- N_tar_vals[which.min(abs(N_tar_vals - N_obs_val))[1]]
}
N_tar_vals_sub[which_repl] <- repl_with

dat_5_pwr <- func_table_pwr(dat)
dat_5_pe <- func_table_pe(dat)

dat_5_pwr
dat_5_pe



# ------------------------------------------------------------------------------
# PROBLEM 6

rm(dat)

fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-02-glmm_testcoef_agg.rds")

dat <- readRDS(fpath_tmp)
head(dat)
table(dat$cnt)
table(dat$name)
table(dat$eff_tar)
sort(unique((dat$N_tar)))

N_tar_vals <- sort(unique(dat$N_tar))
N_obs_val  <- sort(unique(dat$N_obs))[1]
N_tar_idx  <- round(seq(1, length(N_tar_vals), length.out = N_tar_NN))
N_tar_vals_sub <- N_tar_vals[N_tar_idx]
# replace one value with N_obs
which_repl <- which.min(abs(N_tar_vals_sub - N_obs_val))[1]
repl_with <- N_obs_val
if (!(repl_with %in% N_tar_vals)){
  repl_with <- N_tar_vals[which.min(abs(N_tar_vals - N_obs_val))[1]]
}
N_tar_vals_sub[which_repl] <- repl_with

dat_6_pwr <- func_table_pwr(dat)
dat_6_pe <- func_table_pe(dat)

dat_6_pwr
dat_6_pe

library(stargazer)

stargazer(dat_1_pwr, summary = FALSE)
stargazer(dat_2_pwr, summary = FALSE)
stargazer(dat_3_pwr, summary = FALSE)
stargazer(dat_4_pwr, summary = FALSE)
stargazer(dat_5_pwr, summary = FALSE)
stargazer(dat_6_pwr, summary = FALSE)

stargazer(dat_1_pe, summary = FALSE)
stargazer(dat_2_pe, summary = FALSE)
stargazer(dat_3_pe, summary = FALSE)
stargazer(dat_4_pe, summary = FALSE)
stargazer(dat_5_pe, summary = FALSE)
stargazer(dat_6_pe, summary = FALSE)




