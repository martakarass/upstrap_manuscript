#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL/2021-12-03-lm_testcoef_vary_covprop_raw/*
#' 
#' cd $ups/numerical_experiments/R
#' 
#' Rnosave run_lm_testcoef_vary_covprop.R -t 1-1000 -tc 40 -N JOB_lm_testcoef_vary_covprop
#' 
#' ls -l -d *JOB_lm_testcoef_vary_covprop*
#' rm JOB_lm_testcoef_vary_covprop*


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
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-lm_testcoef_vary_covprop_raw")
dir.create(path = res_fdir_raw)

# experiment parameters
# N_obs   <- 50   # size of each of the two "treatment arms" 
N_obs   <- 50
coef_x0 <- 0 
coef_x1 <- 0.5
coef_x2 <- 1
coef_x3 <- -1
# sd_sigma  <- 1
sd_sigma      <- 0.5
N_tar_grid    <- seq(20, 200, by = 30)
N_tar_max     <- max(N_tar_grid)
N_tar_grid_l  <- length(N_tar_grid)
# eff_tar_grid  <- c(0.5)
eff_tar       <- 0.5
eff_tru       <- coef_x1
cov_prop_grid <- c(0.5, 0.3, 0.1, 0.05)

# number of boot repetitions within one experiment, one setup
B_boot  <- 1000
R_powertrue  <- 1000 * 10
# B_boot  <- 3
# R_powertrue  <- 3



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATIONS

mat_out_all <- data.frame()

message(paste0("ESTIMATE POWER WITH upstrap, for fixed target effect size = ", eff_tar))
for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]
  print(paste0("cov_prop = ", cov_prop))
  
  # ------------------------------------------------------------------------------
  # SIMULATE OBSERVED SAMPLE
  
  subjid_arm_i <- rep(1 : N_obs, each = 2)
  subjid_i  <- 1 : (N_obs * 2)   # subject ID unique in data set 
  x1_i      <- rep(c(0, 1), times = N_obs)
  x2_i      <- rbinom(n = N_obs * 2, size = 1, prob = 0.5)
  x3_i      <- runif(n = N_obs * 2, min = 0, max = 1)
  eps_i     <- rnorm(N_obs * 2, sd = sd_sigma)
  y_i       <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) +  (coef_x3 * x3_i) + eps_i
  dat       <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i, 
                          subjid = subjid_i, subjid_arm = subjid_arm_i)
  # indexes of individuals in each two "treatment arms" 
  dat_0_idx <- which(dat$x1 == 0)
  dat_1_idx <- which(dat$x1 == 1)
  # get observed effect 
  fit_obs  <- lm(y ~ x1 + x2 + x3, data = dat)
  
  # ------------------------------------------------------------------------------
  # ESTIMATE POWER WITH upstrap
  
  # update sample to represent target effect (here: no update)
  dat_upd <- dat
  dat_upd$y <- dat_upd$y + (eff_tar - coef(fit_obs)["x1"]) * dat_upd$x1
  # define object to store values across B resamplings
  mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
  ##
  for (bb in 1:B_boot){ # bb <- 100
    for (rr in 1 : N_tar_grid_l){ 
      tryCatch({
        N_tar      <- N_tar_grid[rr]
        size_0 <- round(N_tar * 2 * (1-cov_prop))
        size_1 <- (N_tar * 2) - size_0
        dat_bb_rr_idx_0 <- sample(dat_0_idx, size = size_0, replace = TRUE)
        dat_bb_rr_idx_1 <- sample(dat_1_idx, size = size_1, replace = TRUE)
        dat_bb_rr_idx   <- c(dat_bb_rr_idx_0, dat_bb_rr_idx_1)
        dat_bb_rr       <- dat_upd[dat_bb_rr_idx, ]
        fit_bb_rr  <- lm(y ~ x1 + x2 + x3, data = dat_bb_rr)
        pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
        mat_boot[rr, bb] <- (pval_bb_rr < 0.05) * 1
      }, error = function(e) {message(e)})
    }
  }
  value <- rowMeans(mat_boot, na.rm = TRUE)
  # save to file 
  mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
  mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
  mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
  mat_out_tmp$name          <- "upstrap_power"
  mat_out_tmp$eff_tru       <- eff_tru
  mat_out_tmp$eff_tar       <- eff_tar
  mat_out_tmp$sd_sigma      <- sd_sigma
  mat_out_tmp$value         <- value
  mat_out_tmp$prop_success  <- mean(!is.na(mat_boot))
  mat_out_tmp$cov_prop      <- cov_prop
  mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
  rm(mat_out_tmp, value, mat_boot)
}



# ------------------------------------------------------------------------------
# ADDITIONAL SIMULATIONS TO BE RUN ONCE 

if (arrayjob_idx == 1){
  message(paste0("RUN ADDITIONAL SIMULATIONS TO BE RUN ONCE"))
  
  for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]
    print(paste0("cov_prop = ", cov_prop))
    
    for (bb in 1 : R_powertrue){
      set.seed(bb)
      
      value <- rep(NA, N_tar_grid_l)
      for (rr in 1 : N_tar_grid_l){  
        tryCatch({
          size_0 <- round(N_tar * 2 * (1-cov_prop))
          size_1 <- (N_tar * 2) - size_0  # subjid_arm_i <- rep(1 : N_tar_max, each = 2)
          # subjid_i  <- 1 : (N_tar_max * 2)   # subject ID unique in data set
          x1_i      <- c(rep(0, size_0), rep(1, size_1))
          x2_i      <- rbinom(n = N_tar * 2, size = 1, prob = 0.5)
          x3_i      <- runif(n  = N_tar * 2, min = 0, max = 1)
          eps_i     <- rnorm(n  = N_tar * 2, sd = sd_sigma)
          # use eff_tar
          y_i       <- coef_x0 + (eff_tar * x1_i) + (coef_x2 * x2_i) +  (coef_x3 * x3_i) + eps_i
          dat_bb_rr <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i)
          fit_bb_rr  <- lm(y ~ x1 + x2 + x3, data = dat_bb_rr)
          pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
          value[rr] <- (pval_bb_rr < 0.05) * 1
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
