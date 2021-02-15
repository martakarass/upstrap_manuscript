
#' This script generate
#' - data results
#' showing how sampling distribution of sample mean and sample standard deviation
#' affect the (a) mean, (b) median of power across multiple experiment repetitions. 

rm(list = ls())

library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# dir to save results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-02-10-onesample_ttest_aggegating_comparison")

# simulation parameters
N0     <- 50
N1_min <- N0
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu <- 0.3
sigma2_v <- 1
sigma_v  <- 1
rep_R <- 10000

# factor levels, labels 
aggstat_name_level <- c("power_est_median", "power_est_mean")
aggstat_name_label <- c("median", "mean")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# get theoretical results for sampling distribution of s sample estimator of sd 
# E(s) = \sigma * c_{4}(n)
# sd(s) = \sigma * \sqrt{1-c_{4}^{2}(n)}
s_mean    <- get_mean_s(n = N0, sigma = sigma_v)
s_sd      <- get_sd_s(n = N0, sigma = sigma_v)
xbar_mean <- mu
xbar_sd   <- sigma_v/sqrt(N0)

# simulate draws from sampling distribution of sample mean 
set.seed(1)
xbar_obs_vec   <- rnorm(rep_R, mean = mu, sd = sigma_v/sqrt(N0))
xbar_obs_vec_l <- length(xbar_obs_vec)

# simulate draws from sampling distribution of sample standard deviation (s) 
set.seed(1)
s_obs_vec   <- rnorm(rep_R, mean = s_mean, sd = s_sd)
s_obs_vec_l <- length(s_obs_vec)


# ------------------------------------------------------------------------------
# (1) variate the sample mean estimator 

t1 <- Sys.time()

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_obs_vec[i], sd = s_mean, 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_1 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_1_agg <- 
  plt_df_1 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))


# ------------------------------------------------------------------------------
# (2) variate the sample standard deviation estimator 

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_mean, sd = s_obs_vec[i], 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_2 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_2_agg <- 
  plt_df_2 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))


# ------------------------------------------------------------------------------
# (3) variate the sample mean, (b) sample standard deviation 

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_obs_vec[i], sd = s_obs_vec[i], 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_3 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_3_agg <- 
  plt_df_3 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))

t2 <- Sys.time()
t2 - t1
# 5 min

# sample the idx we store for the purpose of presentation
idx_sub <- seq(1, to = rep_R, length.out = 100)
plt_df_1_sub <- plt_df_1 %>% filter(idx %in% idx_sub)
plt_df_2_sub <- plt_df_2 %>% filter(idx %in% idx_sub)
plt_df_3_sub <- plt_df_3 %>% filter(idx %in% idx_sub)
saveRDS(plt_df_1_sub, file = paste0(results_dir, "/plt_df_1_sub.rds"))
saveRDS(plt_df_2_sub, file = paste0(results_dir, "/plt_df_2_sub.rds"))
saveRDS(plt_df_3_sub, file = paste0(results_dir, "/plt_df_3_sub.rds"))
saveRDS(plt_df_1_agg, file = paste0(results_dir, "/plt_df_1_agg.rds"))
saveRDS(plt_df_2_agg, file = paste0(results_dir, "/plt_df_2_agg.rds"))
saveRDS(plt_df_3_agg, file = paste0(results_dir, "/plt_df_3_agg.rds"))
