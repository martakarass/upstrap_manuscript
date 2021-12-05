#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL/2021-12-03-twosample_ttest_vary_covprop_raw/*
#' 
#' cd $ups/numerical_experiments/R
#' 
#' Rnosave run_twosample_ttest_vary_covprop.R -t 1-1000 -tc 40 -N JOB_twosample_ttest_vary_covprop
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
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-twosample_ttest_vary_covprop_raw")
dir.create(path = res_fdir_raw)

# experiment parameters
N_obs         <- 50 # size of each of the two "treatment arms" 
N_tar_grid    <- seq(20, 200, by = 30)
N_tar_max     <- max(N_tar_grid)
N_tar_grid_l  <- length(N_tar_grid)
sd_sigma      <- 0.5
eff_tar       <- 0.3
eff_tru       <- 0.3
cov_prop_grid <- c(0.5, 0.3, 0.1, 0.05)

# number of boot repetitions within one experiment, one setup
# B_boot  <- 1000
# R_powertrue  <- 1000 * 10
B_boot  <- 3
R_powertrue  <- 3


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATIONS

message(paste0("ESTIMATE POWER WITH upstrap, for fixed target effect size = ", eff_tar))

mat_out_all <- data.frame()
for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]
  print(paste0("cov_prop = ", cov_prop))

  # ------------------------------------------------------------------------------
  # SIMULATE OBSERVED SAMPLE
  sample_1_i <- rnorm(n = N_obs, mean = 0,       sd = sd_sigma)
  sample_2_i <- rnorm(n = N_obs, mean = eff_tru, sd = sd_sigma)

  # ------------------------------------------------------------------------------
  # ESTIMATE POWER WITH upstrap
  
  # update the outcome in the sample to represent the target effect size 
  sample_i_meandiff   <-  mean(sample_2_i) - mean(sample_1_i)
  sample_1_i_upd <- sample_1_i 
  sample_2_i_upd <- sample_2_i + (eff_tar - sample_i_meandiff) * 1
  # define object to store values across B resamplings
  mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
  for (bb in 1:B_boot){ 
    for (rr in 1 : N_tar_grid_l){ 
      tryCatch({
        N_tar          <- N_tar_grid[rr]
        sample_1_bb_rr <- sample(sample_1_i_upd, size = round(N_tar * 2 * (1-cov_prop)),     replace = TRUE)
        sample_2_bb_rr <- sample(sample_2_i_upd, size = round(N_tar * 2 * cov_prop), replace = TRUE)
        mat_boot[rr, bb] <- (t.test(x = sample_1_bb_rr, y = sample_2_bb_rr, var.equal = TRUE)$p.value < 0.05)
      }, error = function(e) {message(e)})
    }
  }
  value <- rowMeans(mat_boot, na.rm = TRUE)
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
  mat_out_tmp$cov_prop      <- cov_prop
  mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
  rm(mat_out_tmp, value)
}


# ------------------------------------------------------------------------------
# ADDITIONAL SIMULATIONS TO BE RUN ONCE 

if (arrayjob_idx == 1){
  message(paste0("RUN ADDITIONAL SIMULATIONS TO BE RUN ONCE"))
  
  for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]
    print(paste0("cov_prop = ", cov_prop))
    for (bb in 1 : R_powertrue){
      set.seed(bb)
      sample_1_i <- rnorm(n = round(N_tar_max * 2 * (1-cov_prop)),  mean = 0,       sd = sd_sigma)
      sample_2_i <- rnorm(n = round(N_tar_max * 2 * cov_prop),      mean = eff_tru, sd = sd_sigma)
      value <- rep(NA, N_tar_grid_l)
      for (rr in 1 : N_tar_grid_l){  # rr <- 10
        tryCatch({
          N_tar         <- N_tar_grid[rr]
          sample_1_bb_rr <- sample_1_i[1 : round(N_tar * 2 * (1-cov_prop))]
          sample_2_bb_rr <- sample_2_i[1 : round(N_tar * 2 * cov_prop)]
          value[rr] <- (t.test(x = sample_1_bb_rr, y = sample_2_bb_rr, var.equal = TRUE)$p.value < 0.05)
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
      mat_out_tmp$cov_prop      <- cov_prop
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
