
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)
library(stargazer)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

my_pal <- pal_futurama()(12)
show_col(my_pal)
show_col(my_pal[c(4,9)])


# TABLE 1a: one-sample t-test power  ---------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-onesample_ttest_agg/")
onesample_power_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_power_est_bootCI.rds"))
onesample_samplesize_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_samplesize_est_bootCI.rds"))
onesample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret.rds")) %>%
  filter(power_est >= 0.8) %>% filter(N1 == min(N1)) 
  

tbl_1a <- 
  onesample_samplesize_df %>%
  mutate(val = paste0(power_est_aggmedian, " [", power_est_aggmedian_lwr, ", ", 
                      round(power_est_aggmedian_upr, 1), "]")) %>%
  select(method_name, N0, val) %>%
  pivot_wider(names_from = method_name, values_from = val) %>%
  mutate(N1_true = rep(onesample_power_theoret_df$N1, 3), .after = N0) %>%
  mutate(problem = rep("one-sample t-test", 3), .before = N0)


# TABLE 1b: two-sample t-test power  ---------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-twosample_ttest_agg")
twosample_power_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_power_est_bootCI.rds"))
twosample_samplesize_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_samplesize_est_bootCI.rds"))
twosample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret")) %>%
  filter(power_est >= 0.8) %>% filter(N1 == min(N1)) 

tbl_1b <- 
  twosample_samplesize_df %>%
  mutate(val = paste0(power_est_aggmedian, " [", power_est_aggmedian_lwr, ", ", 
                      round(power_est_aggmedian_upr, 1), "]")) %>%
  select(method_name, N0, val) %>%
  pivot_wider(names_from = method_name, values_from = val) %>%
  arrange(N0) %>%
  mutate(N1_true = rep(twosample_power_theoret_df$N1, 3), .after = N0) %>%
  mutate(problem = rep("two-sample t-test", 3), .before = N0)

tbl_out <- 
  rbind(tbl_1a, tbl_1b) %>% 
  as.data.frame() %>%
  rename(N_observed = N0, M_theoret = N1_true, 
         M_est_powerttest = powerttest, M_est_upstrap = upstrap)

stargazer::stargazer(tbl_out, summary = FALSE)



