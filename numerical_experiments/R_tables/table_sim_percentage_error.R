
#' Compute table summary of PE 

rm(list = ls())
library(tidyverse)
library(here)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# READ AND PROCESS DATA 

N_tar_NN <- 6


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

dat_1_pwr <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "upstrap_power") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 1, .before = everything())
dat_1_pwr

dat_1_pe <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "ups_cmp_power_pe") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 1, .before = everything())
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

dat_2_pwr <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "upstrap_power") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
dat_2_pwr

dat_2_pe <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "ups_cmp_power_pe") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
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

dat_3_pwr <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "upstrap_power") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
dat_3_pwr

dat_3_pe <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "ups_cmp_power_pe") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
dat_3_pe



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

dat_5_pwr <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "upstrap_power") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
dat_5_pwr

dat_5_pe <- 
  dat %>%
  filter(N_tar %in% N_tar_vals_sub) %>%
  filter(name == "ups_cmp_power_pe") %>%
  mutate(value_mean_f   = sprintf("%.2f", value_mean)) %>%
  mutate(value_sd_f     = sprintf("%.2f", value_sd)) %>%
  mutate(value_f        = paste0(value_mean_f, " (", value_sd_f, ")")) %>%
  select(N_obs, N_tar, eff_tar, value_f) %>%
  pivot_wider(names_from = eff_tar, values_from = value_f) %>%
  mutate(problem_id = 2, .before = everything())
dat_5_pe





