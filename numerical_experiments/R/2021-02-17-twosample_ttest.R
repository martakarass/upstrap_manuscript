
#' This script estimate power of rejecting H0 in two-sample t-test problem. 

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_agg <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-twosample_ttest_agg")
res_fdir_raw <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-twosample_ttest_raw")
# create dirs if any does not exist
dir.create(path = res_fdir_agg)
dir.create(path = res_fdir_raw)
message(paste0("dir.exists(path = res_fdir_agg): ", dir.exists(path = res_fdir_agg)))
message(paste0("dir.exists(path = res_fdir_raw): ", dir.exists(path = res_fdir_raw)))

# x   y
# 1 0.10  12
# 2 0.30  47
# 3 0.50  87
# 4 0.70 139
# 5 0.90 235
# 6 0.95 290

# experiment parameteres
N0_grid <- c(47, 87, 139)
N1_min <- 5
N1_max <- 300 
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
out <- sapply(N1_grid, function(n_tmp) power.t.test(n = n_tmp, delta = mu, sd = 1, type = "two.sample")$power)
out_df_1 <- data.frame(N1 = N1_grid, power_est = out)
out_df_fpath <- paste0(res_fdir_agg, "/res_theoret")
saveRDS(out_df_1, out_df_fpath)


# compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
cum_rejectH0_twosample_ttest <- function(vals1, vals2){
  vals_cumnx      <- seq_along(vals1)
  vals_cumdf      <- 2 * vals_cumnx - 2
  vals1_cummean    <- cumsum(vals1) / vals_cumnx
  vals2_cummean    <- cumsum(vals2) / vals_cumnx
  vals_cummeandiff <- vals1_cummean - vals2_cummean
  vals1_sumsquares <- (cumsum(vals1 ^ 2) - cumsum(vals1) ^ 2 / seq_along(vals1))
  vals2_sumsquares <- (cumsum(vals2 ^ 2) - cumsum(vals2) ^ 2 / seq_along(vals2))
  vals_cumvarpool  <- (vals1_sumsquares + vals2_sumsquares) / vals_cumdf
  vals_cumsizeinvsum  <- (1 / vals_cumnx) + (1 / vals_cumnx)
  vals_cumstderr   <- sqrt(vals_cumvarpool * vals_cumsizeinvsum)
  vals_cumststat   <- vals_cummeandiff / vals_cumstderr
  vals_cumpval     <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
  vals_reject_H0   <- vals_cumpval < 0.05
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
  x_mat_gr1  <- matrix(rnorm(n = rep_n * N0, mean = 0), nrow = rep_n)
  x_mat_gr2  <- matrix(rnorm(n = rep_n * N0, mean = mu), nrow = rep_n)
  mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  for (i in 1:rep_n){
    sample_i_gr1 <- x_mat_gr1[i, ]
    sample_i_gr2 <- x_mat_gr2[i, ]
    # sample_i: generate power estimate with power.t.test()
    sample_i_meandiff   <-  mean(sample_i_gr1) - mean(sample_i_gr2)
    sample_i_var_pooled <- (var(sample_i_gr1) * (N0 - 1) + var(sample_i_gr2) * (N0 - 1)) / (N0 + N0 - 2)
    sample_i_sd_pooled  <- sqrt(sample_i_var_pooled)
    mat_out_powerttest[i, ] <- sapply(N1_grid, function(N1_tmp){
      (power.t.test(n = N1_tmp, 
                    delta = sample_i_meandiff, 
                    sd = sample_i_sd_pooled, 
                    sig.level = 0.05, type = "two.sample", alternative = "two.sided"))$power
    })
    # sample_i: generate power estimate with upstrap()
    boot_resamples_i_gr1 <- matrix(sample(x = sample_i_gr1, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_gr2 <- matrix(sample(x = sample_i_gr2, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_rejectH0 <- lapply(1:B_boot, function(k) cum_rejectH0_twosample_ttest(
      boot_resamples_i_gr1[k, ],
      boot_resamples_i_gr2[k, ]
    ))
    boot_resamples_i_rejectH0 <- do.call(rbind, boot_resamples_i_rejectH0)
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



