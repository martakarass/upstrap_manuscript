
#' This script estimate power of rejecting H0 in one-sample t-test problem. 
#' 
#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL/2021-08-07-twosample_ttest_raw/*
#' cd $ups/numerical_experiments/R
#' Rnosave run_twosample_ttest.R -N JOB_twosample
#' 
#' ls -l -d *JOB_twosample*
#' rm JOB_twosample*


rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-twosample_ttest_raw")
dir.create(path = res_fdir_raw)

# experiment parameters
N_obs_grid    <- c(50, 100)
N_tar_grid    <- seq(20, 200, by = 30)
N_tar_max     <- max(N_tar_grid)
N_tar_grid_l  <- length(N_tar_grid)
eff_tar_grid  <- c(0.3, 0.4)

# data generating model
eff_tru <- 0.3
simga2 <- 1

# number of repetitions of experiment 
R_rep   <- 1000
# number of boot repetitions within one experiment, one setup
B_boot  <- 1000


# ------------------------------------------------------------------------------
# HELPER FUNCTIONS

# function to compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
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


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATIONS

for (arrayjob_idx in 1 : R_rep){ # arrayjob_idx <- 1
  set.seed(arrayjob_idx)
  message(paste0("arrayjob_idx: ", arrayjob_idx))
  mat_out_all <- data.frame()
  
  # iterate over observed sample size N_obs
  for (N_obs in N_obs_grid){ # N_obs <- N_obs_grid[1]
  
    # simulate sample observed in this experiment repetition
    sample_1_i <- rnorm(n = N_obs, mean = 0)
    sample_2_i <- rnorm(n = N_obs, mean = eff_tru)
    # observed quantities
    sample_i_meandiff   <-  mean(sample_2_i) - mean(sample_1_i)
    sample_i_var_pooled <- (var(sample_1_i) * (N_obs - 1) + var(sample_2_i) * (N_obs - 1)) / (N_obs + N_obs - 2)
    sample_i_sd_pooled  <- sqrt(sample_i_var_pooled)
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH power.t.test(), for observed effect size 
    
    value <- sapply(N_tar_grid, function(n_tmp){
      out_test <- power.t.test(n = n_tmp,  delta = sample_i_meandiff,  sd = sample_i_sd_pooled, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
      out_test$power
    })
    mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
    mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
    mat_out_tmp$name          <- "powerttest_power"
    mat_out_tmp$eff_tru       <- eff_tru
    mat_out_tmp$eff_tar       <- NA # observed power case
    mat_out_tmp$value         <- value
    mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
    
    
    # ESTIMATE POWER WITH power.t.test(), for predefined effect size 
    for (eff_tar in eff_tar_grid){ # eff_tar <- eff_tar_grid[1]
      value <- sapply(N_tar_grid, function(n_tmp){
        out_test <- power.t.test(n = n_tmp,  delta = eff_tar,  sd = sample_i_sd_pooled, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
        out_test$power
      })
      mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
      mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
      mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
      mat_out_tmp$name          <- "powerttest_power"
      mat_out_tmp$eff_tru       <- eff_tru
      mat_out_tmp$eff_tar       <- eff_tar
      mat_out_tmp$value         <- value
      mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
      rm(mat_out_tmp, value)
    }
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH upstrap, for observed effect size 
    
    sample_1_i_upd <- sample_1_i 
    sample_2_i_upd <- sample_2_i # no update 
    boot_resamples_1_i <- matrix(sample(x = sample_1_i_upd, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
    boot_resamples_2_i <- matrix(sample(x = sample_2_i_upd, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
    boot_resamples_i_rejectH0 <- lapply(1 : B_boot, function(k) cum_rejectH0_twosample_ttest(
      boot_resamples_1_i[k, ],
      boot_resamples_2_i[k, ]
    ))
    boot_resamples_i_rejectH0 <- do.call(rbind, boot_resamples_i_rejectH0)
    value <- apply(boot_resamples_i_rejectH0[, N_tar_grid], 2, mean, na.rm = TRUE)
    ##
    mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
    mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
    mat_out_tmp$name          <- "upstrap_power"
    mat_out_tmp$eff_tru       <- eff_tru
    mat_out_tmp$eff_tar       <- NA
    mat_out_tmp$value         <- value
    mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
    
    
    # ESTIMATE POWER WITH upstrap, for predefined effect size 
    for (eff_tar in eff_tar_grid){
      # update the outcome in the sample to represent the target effect size 
      sample_1_i_upd <- sample_1_i 
      sample_2_i_upd <- sample_2_i + (eff_tar - sample_i_meandiff) * 1
      boot_resamples_1_i <- matrix(sample(x = sample_1_i_upd, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
      boot_resamples_2_i <- matrix(sample(x = sample_2_i_upd, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
      boot_resamples_i_rejectH0 <- lapply(1 : B_boot, function(k) cum_rejectH0_twosample_ttest(
        boot_resamples_1_i[k, ],
        boot_resamples_2_i[k, ]
      ))
      boot_resamples_i_rejectH0 <- do.call(rbind, boot_resamples_i_rejectH0)
      value <- apply(boot_resamples_i_rejectH0[, N_tar_grid], 2, mean, na.rm = TRUE)
      ##
      mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
      mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
      mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
      mat_out_tmp$name          <- "upstrap_power"
      mat_out_tmp$eff_tru       <- eff_tru
      mat_out_tmp$eff_tar       <- eff_tar
      mat_out_tmp$value         <- value
      mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
      rm(mat_out_tmp, value)
    }
  }
  
  # ------------------------------------------------------------------------------
  # ESTIMATE "TRUE POWER"
  
  for (eff_tar in eff_tar_grid){
    value <- numeric(N_tar_grid_l)
    for (rr in 1 : N_tar_grid_l){
      N_obs <- N_tar_grid[rr]
      sample_1_i <- rnorm(n = N_obs, mean = 0)
      sample_2_i <- rnorm(n = N_obs, mean = eff_tar)
      value[rr] <- (t.test(x = sample_1_i, y = sample_2_i, var.equal = TRUE)$p.value < 0.05)
    }
    mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
    mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
    mat_out_tmp$name          <- "true_power"
    mat_out_tmp$eff_tru       <- NA
    mat_out_tmp$eff_tar       <- eff_tar
    mat_out_tmp$value         <- value
    mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
  }
  
  
  # ------------------------------------------------------------------------------
  # SAVE TO FILE 
  out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
  saveRDS(object = mat_out_all, file = out_fpath_raw)
}