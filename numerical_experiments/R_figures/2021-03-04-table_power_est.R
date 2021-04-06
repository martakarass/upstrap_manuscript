
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)
library(stargazer)

source(here::here("numerical_experiments/R/config_utils.R"))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TABLE 1a: one-sample t-test power  -------------------------------------------

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


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TABLE 2a: LM    --------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-03-01-lm_trt_agg/")
samplesize_df <- readRDS(paste0(fdir_tmp, "/lm_trt_samplesize_est_bootCI.rds"))

tbl_2a_gold <- 
  samplesize_df %>% filter(method_name == "out_gold_LM") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                 round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2a_upst <- 
  samplesize_df %>% filter(method_name == "out_boot_LM") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                      round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2a <- 
  tbl_2a_gold %>% 
  left_join(tbl_2a_upst, by = "N0") %>%
  mutate(problem = rep("Problem (3)", 3), .before = everything())
tbl_2a


# TABLE 2b: GLM    -------------------------------------------------------------
fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-03-01-glm_trt_agg/")
samplesize_df <- readRDS(paste0(fdir_tmp, "/glm_trt_samplesize_est_bootCI.rds"))

tbl_2b_gold <- 
  samplesize_df %>% filter(method_name == "out_gold_GLM") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                      round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2b_upst <- 
  samplesize_df %>% filter(method_name == "out_boot_GLM") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                     round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2b <- 
  tbl_2b_gold %>% 
  left_join(tbl_2b_upst, by = "N0") %>%
  mutate(problem = rep("Problem (4)", 3), .before = everything())
tbl_2b


# TABLE 2c: LMMM    ------------------------------------------------------------
fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-lmm_trt_agg/")
samplesize_df <- readRDS(paste0(fdir_tmp, "/lmm_trt_samplesize_est_bootCI.rds"))

tbl_2c_gold_GEE <- 
  samplesize_df %>% filter(method_name == "out_gold_GEE") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                      round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2c_boot_GEE <- 
  samplesize_df %>% filter(method_name == "out_boot_GEE") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                     round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2c_GEE <- 
  tbl_2c_gold_GEE %>% 
  left_join(tbl_2c_boot_GEE, by = "N0") %>%
  mutate(problem = "Problem (5) -- GEE", .before = everything())
tbl_2c_GEE


tbl_2c_gold_LMM <- 
  samplesize_df %>% filter(method_name == "out_gold_LMM") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                      round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2c_boot_LMM <- 
  samplesize_df %>% filter(method_name == "out_boot_LMM") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                     round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2c_LMM <- 
  tbl_2c_gold_LMM %>% 
  left_join(tbl_2c_boot_LMM, by = "N0") %>%
  mutate(problem = "Problem (5) -- LMM", .before = everything())
tbl_2c_LMM

tbl_2c <- rbind(tbl_2c_GEE, tbl_2c_LMM)
tbl_2c


# TABLE 2d: LMMM    ------------------------------------------------------------
fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-glmm_trt_agg/")
samplesize_df <- readRDS(paste0(fdir_tmp, "/glmm_trt_samplesize_est_bootCI.rds"))

tbl_2d_gold_GEE <- 
  samplesize_df %>% filter(method_name == "out_gold_GEE") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                      round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2d_boot_GEE <- 
  samplesize_df %>% filter(method_name == "out_boot_GEE") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                     round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2d_GEE <- 
  tbl_2d_gold_GEE %>% 
  left_join(tbl_2d_boot_GEE, by = "N0") %>%
  mutate(problem = "Problem (6) -- GEE", .before = everything())
tbl_2d_GEE


tbl_2d_gold_GLMM <- 
  samplesize_df %>% filter(method_name == "out_gold_GLMM") %>%
  mutate(
    out_gold = paste0(est_mean, " [", round(est_bootci_mean_lwr, 1), ", ", 
                      round(est_bootci_mean_upr, 1), "]")) %>%
  select(N0, out_gold) 
tbl_2d_boot_GLMM <- 
  samplesize_df %>% filter(method_name == "out_boot_GLMM") %>%
  mutate(
    out_ups = paste0(est_median, " [", round(est_bootci_median_lwr, 1), ", ", 
                     round(est_bootci_median_upr, 1), "]")) %>%
  select(N0, out_ups)
tbl_2d_GLMM <- 
  tbl_2d_gold_GLMM %>% 
  left_join(tbl_2d_boot_GLMM, by = "N0") %>%
  mutate(problem = "Problem (6) -- GLMM", .before = everything())
tbl_2d_GLMM

tbl_2d <- rbind(tbl_2d_GEE, tbl_2d_GLMM)
tbl_2d


# TABLE 2: FINAL   -------------------------------------------------------------

tbl_2 <- rbind(
  tbl_2a,
  tbl_2b,
  tbl_2c,
  tbl_2d
)
tbl_2

stargazer::stargazer(tbl_2, summary = FALSE)




