
#' This script estimate power of rejecting H0 in one-sample t-test problem. 

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-15-onesample_ttest_agg")
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-15-onesample_ttest_raw")
# create dirs if any does not exist
dir.create(path = res_fdir_agg)
dir.create(path = res_fdir_raw)
message(paste0("dir.exists(path = res_fdir_agg): ", dir.exists(path = res_fdir_agg)))
message(paste0("dir.exists(path = res_fdir_raw): ", dir.exists(path = res_fdir_raw)))

# experiment parameters
N0_grid <- c(30, 50, 100, 150)
N1_min <- 30
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu     <- 0.3
simga2 <- 1

# number of repetitions of experiment 
rep_n   <- 10000
# number of boot repetitions within one experiment, one setup
B_boot  <- 1000


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# get the estimates using theoretical (gold standard) results
out <- sapply(N1_grid, function(n_tmp) power.t.test(n = n_tmp, delta = mu, sd = 1, type = "one.sample")$power)
out_df_1 <- data.frame(N1 = N1_grid, power_est = out)
out_df_fpath <- paste0(res_fdir_agg, "/res_theoret.rds")
saveRDS(out_df_1, out_df_fpath)

# function to compute cumulative var
cumvar <- function (x, sd = FALSE) {
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  # if (sd) v <- sqrt(v)
  v
}

# function to compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
vals_cum_reject_H0 <- function(vals){
  vals_cumnx      <- seq_along(vals)
  vals_cumdf      <- vals_cumnx - 1
  vals_cummean    <- cumsum(vals) / vals_cumnx
  vals_cumvar     <- cumvar(vals)
  vals_cumstderr  <- sqrt(vals_cumvar / vals_cumnx)
  vals_cumststat  <- vals_cummean / vals_cumstderr
  vals_cumpval    <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
  vals_reject_H0  <- vals_cumpval < 0.05
  return(vals_reject_H0)
}

N0_vec        <- numeric()
N1_vec        <- numeric()
group_vec     <- numeric()
power_est_vec <- numeric()

t1 <- Sys.time()
set.seed(123)
for (N0 in N0_grid){ # N0 <- 50; i <- 1
  message(N0)
  N1_grid <- N0 : N1_max
  N1_grid_l <- length(N1_grid)
  # set of samples for current N0
  x_mat <- matrix(rnorm(n = rep_n * N0, mean = mu), nrow = rep_n)
  mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  for (i in 1:rep_n){
    sample_i <- x_mat[i, ]
    # sample_i: generate power estimate with power.t.test()
    sample_i_mean <- mean(sample_i)
    sample_i_sd   <- sd(sample_i)
    mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
      (power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
    })
    # sample_i: generate power estimate with upstrap()
    boot_resamples_i <- matrix(sample(x = sample_i, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
    mat_out_upstrap[i, ] <- apply(boot_resamples_i_rejectH0[, N1_grid], 2, mean, na.rm = TRUE)
  }
  
  # aggregate raw data
  out_powerttest_mean   <- apply(mat_out_powerttest, 2, mean)
  out_powerttest_median <- apply(mat_out_powerttest, 2, median)
  out_upstrap_mean      <- apply(mat_out_upstrap, 2, mean)
  out_upstrap_median    <- apply(mat_out_upstrap, 2, median)
  
  # store raw data aggregates 
  N0_vec        <- c(N0_vec, rep(N0, 4 * N1_grid_l))
  N1_vec        <- c(N1_vec, rep(N1_grid, 4))
  group_vec     <- c(group_vec, 
                     rep("powerttest_mean", N1_grid_l), 
                     rep("powerttest_median", N1_grid_l),
                     rep("upstrap_mean", N1_grid_l), 
                     rep("upstrap_median", N1_grid_l))
  power_est_vec <- c(power_est_vec, 
                     out_powerttest_mean,
                     out_powerttest_median,
                     out_upstrap_mean,
                     out_upstrap_median)
  
  # save raw data 
  saveRDS(mat_out_powerttest, paste0(res_fdir_raw, "/mat_out_powerttest_N0_", N0, ".rds"))
  saveRDS(mat_out_upstrap,    paste0(res_fdir_raw, "/mat_out_upstrap_N0_", N0, ".rds"))
}

# save raw data aggregates 
out_df_2 <- data.frame(N0 = N0_vec, 
                       N1 = N1_vec, 
                       group = group_vec, 
                       power_est = power_est_vec) 
out_df_fpath <- paste0(res_fdir_agg, "/res_repn_", rep_n, "_boot_", B_boot, "_powerttest_upstrap.rds")
saveRDS(out_df_2, out_df_fpath)

t2 <- Sys.time()
message(t2 - t1)

