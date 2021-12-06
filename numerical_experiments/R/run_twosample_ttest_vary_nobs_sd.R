#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL/2021-12-03-twosample_ttest_vary_nobs_sd_raw/*
#' 
#' cd $ups/numerical_experiments/R
#' 
#' Rnosave run_twosample_ttest_vary_nobs_sd.R -t 1-1000 -tc 40 -N JOB_twosample_ttest_vary_nobs_sd
#' 
#' ls -l -d *JOB_twosample*
#' rm JOB_twosample*


arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# rm(list = ls()); arrayjob_idx <- 1
set.seed(arrayjob_idx)
message(paste0("arrayjob_idx: ", arrayjob_idx))

library(here)
library(tidyverse)
library(matrixStats)
library(simr)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-twosample_ttest_vary_nobs_sd_raw")
dir.create(path = res_fdir_raw)

# experiment parameters
# N_obs   <- 50   # size of each of the two "treatment arms" 
N_obs_grid    <- c(5, 10, 20, 50)
N_tar_grid    <- seq(20, 200, by = 30)
N_tar_max     <- max(N_tar_grid)
N_tar_grid_l  <- length(N_tar_grid)
sd_sigma_grid <- c(0.4, 0.7, 1)
eff_tar       <- 0.3
eff_tru       <- 0.3

# number of boot repetitions within one experiment, one setup
B_boot  <- 1000
R_powertrue  <- 1000 * 10
# B_boot  <- 3
# R_powertrue  <- 3

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

mat_out_all <- data.frame()

message(paste0("ESTIMATE POWER WITH upstrap, for fixed target effect size = ", eff_tar))
for (N_obs in N_obs_grid){
  for (sd_sigma in sd_sigma_grid){
    print(paste0("N_obs = ", N_obs, ", sd_sigma = ", sd_sigma))
    # N_obs  <- N_obs_grid[1]; sd_sigma <- sd_sigma_grid[2]
    
    # ------------------------------------------------------------------------------
    # SIMULATE OBSERVED SAMPLE
    sample_1_i <- rnorm(n = N_obs, mean = 0, sd = sd_sigma)
    sample_2_i <- rnorm(n = N_obs, mean = eff_tru, sd = sd_sigma)
    # observed quantities
    sample_i_meandiff   <-  mean(sample_2_i) - mean(sample_1_i)
    sample_i_var_pooled <- (var(sample_1_i) * (N_obs - 1) + var(sample_2_i) * (N_obs - 1)) / (N_obs + N_obs - 2)
    sample_i_sd_pooled  <- sqrt(sample_i_var_pooled)
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH upstrap
    
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
    mat_out_tmp$eff_tar       <- eff_tar
    mat_out_tmp$sd_sigma      <- sd_sigma
    mat_out_tmp$value         <- value
    mat_out_tmp$prop_success  <- as.numeric(!is.na(value))
    mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
  }
}


# ------------------------------------------------------------------------------
# ADDITIONAL SIMULATIONS TO BE RUN ONCE 

if (arrayjob_idx == 1){
  message(paste0("RUN ADDITIONAL SIMULATIONS TO BE RUN ONCE"))
  for (sd_sigma in sd_sigma_grid){
    print(paste0("sd_sigma = ", sd_sigma))
    for (bb in 1 : R_powertrue){
      set.seed(bb)
      sample_1_i <- rnorm(n = N_tar_max, mean = 0,       sd = sd_sigma)
      sample_2_i <- rnorm(n = N_tar_max, mean = eff_tru, sd = sd_sigma)
      value <- rep(NA, N_tar_grid_l)
      for (rr in 1 : N_tar_grid_l){  # rr <- 10
        tryCatch({
          N_tar         <- N_tar_grid[rr]
          sample_1_i_rr <- sample_1_i[1 : N_tar]
          sample_2_i_rr <- sample_2_i[1 : N_tar]
          value[rr] <- (t.test(x = sample_1_i_rr, y = sample_2_i_rr, var.equal = TRUE)$p.value < 0.05)
        }, error = function(e) {message(e)})
      }
      mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
      mat_out_tmp$N_obs         <- N_tar_grid
      mat_out_tmp$arrayjob_idx  <- rep(bb, N_tar_grid_l)
      mat_out_tmp$name          <- "true_power"
      mat_out_tmp$eff_tru       <- eff_tar
      mat_out_tmp$eff_tar       <- eff_tar
      mat_out_tmp$sd_sigma      <- sd_sigma
      mat_out_tmp$value         <- value
      mat_out_tmp$prop_success  <- as.numeric(!is.na(value))
      mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
      rm(mat_out_tmp, value)
    }
  }
}



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SAVE TO FILE 

message(paste0("*** COMPLETED ***"))

out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
saveRDS(object = mat_out_all, file = out_fpath_raw)
