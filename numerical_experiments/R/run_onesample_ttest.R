
#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL/2021-08-07-onesample_ttest_raw/*
#' cd numerical_experiments/R
#' Rnosave run_onesample_ttest.R -N JOB_onesample
#' 
#' ls -l -d *JOB_onesample*
#' rm JOB_onesample*

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-onesample_ttest_raw")
# create dirs if any does not exist
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
    sample_i <- rnorm(n = N_obs, mean = eff_tru, sd = sqrt(simga2))
    sample_i_mean <- mean(sample_i)
    sample_i_sd   <- sd(sample_i)
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH power.t.test(), for observed effect size 
    
    value <- sapply(N_tar_grid, function(n_tmp){
      out_test <- power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
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
    for (eff_tar in eff_tar_grid){
      value <- sapply(N_tar_grid, function(n_tmp){
        out_test <- power.t.test(n = n_tmp, delta = eff_tar, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
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
    
    sample_i_updated <- sample_i 
    boot_resamples_i <- matrix(sample(x = sample_i_updated, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
    boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
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
      sample_i_updated <- sample_i + (eff_tar - sample_i_mean)
      boot_resamples_i <- matrix(sample(x = sample_i_updated, size = (B_boot * N_tar_max), replace = TRUE), nrow = B_boot, ncol = N_tar_max)
      boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
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
      sample_i <- rnorm(n = N_obs, mean = eff_tar, sd = sqrt(simga2))
      value[rr] <- (t.test(sample_i)$p.value < 0.05)
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

