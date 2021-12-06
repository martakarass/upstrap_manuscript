#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' cd numerical_experiments/R
#' rm $ups/numerical_experiments/results_CL/2021-12-03-lmm_testcoef_vary_covprop_raw/*
#' Rnosave run_lmm_testcoef_vary_covprop.R -t 1-1000 -tc 75 -N JOB_lmm_testcoef_vary_covprop
#' 
#' ls -l -d *JOB_lmm_testcoef_vary_covprop*
#' rm JOB_lmm_testcoef_vary_covprop*


arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# rm(list = ls()); arrayjob_idx <- 1
set.seed(arrayjob_idx)
message(paste0("arrayjob_idx: ", arrayjob_idx))

library(here)
library(tidyverse)
library(matrixStats)
library(simr)
library(lme4)
library(lmerTest)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-lmm_testcoef_vary_covprop_raw")
# create dirs if any does not exist
dir.create(path = res_fdir_raw)

# experiment parameters
# N_obs   <- 50 # number of subjects in each of the two "treatment arms" 
N_obs   <- 50
ni      <- 3  # number of observations per subject 
coef_x0 <- 0 
coef_x1 <- 0.5
coef_x2 <- 1
coef_x3 <- -1
# sd_tau    <- 1
sd_sigma  <- 0.5
sd_tau    <- 0.5
N_tar_grid    <- seq(20, 200, by = 30)
N_tar_max     <- max(N_tar_grid)
N_tar_grid_l  <- length(N_tar_grid)
# eff_tar_grid  <- c(0.5, 1)
eff_tar  <- 0.5
eff_tru  <- coef_x1
cov_prop_grid <- c(0.5, 0.3, 0.1, 0.05)

# number of boot repetitions within one experiment, one setup
B_boot  <- 1000
R_powertrue <- 1000 * 10
# B_boot  <- 3
# R_powertrue <- 3


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATIONS

# N_obs <- 50; sd_sigma <- 0.5; sd_tau <- sd_sigma

mat_out_all <- data.frame()

for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]
  print(paste0("cov_prop = ", cov_prop))
  
  # ------------------------------------------------------------------------------
  # SIMULATE OBSERVED SAMPLE
  
  # subject-specific
  x1_i      <- rep(c(0, 1), times = N_obs)
  x2_i      <- rbinom(n = N_obs * 2, size = 1, prob = 0.5)
  x3_i      <- runif(n = N_obs * 2, min = 0, max = 1)
  b0_i      <- rnorm(n = (N_obs * 2), mean = 0, sd = sd_tau)
  subjid_i     <- 1 : (N_obs * 2)           # subject ID unique in data set 
  subjid_arm_i <- rep(1 : N_obs, each = 2)  # subject ID unique within "treatment arm" 
  # observation-specific
  x1_ij     <- rep(x1_i, each = ni) 
  x2_ij     <- rep(x2_i, each = ni) 
  x3_ij     <- rep(x3_i, each = ni) 
  b0_ij     <- rep(b0_i, each = ni) 
  eps_ij    <- rnorm(n = (N_obs * 2 * ni), mean = 0, sd = sd_sigma)
  subjid_ij     <- rep(subjid_i, each = ni)            
  subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
  rm(x1_i, x2_i, x3_i, b0_i, subjid_i, subjid_arm_i)
  # data set 
  y_ij      <- (b0_ij + coef_x0) + (coef_x1 * x1_ij) + (coef_x2 * x2_ij) +  (coef_x3 * x3_ij) + eps_ij
  dat       <- data.frame(y = y_ij, x1 = x1_ij, x2 = x2_ij,  x3 = x3_ij,
                          subjid = subjid_ij, subjid_arm = subjid_arm_ij)
  # indexes of individuals in each two "treatment arms" 
  dat_subjid_x1_is0 <- unique(dat %>% filter(x1 == 0) %>% pull(subjid))
  dat_subjid_x1_is1 <- unique(dat %>% filter(x1 == 1) %>% pull(subjid))
  # get observed effect 
  fit_obs  <- lmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat) 
  
  
  # ------------------------------------------------------------------------------
  # ESTIMATE POWER WITH upstrap, for fixed target effect size
  message(paste0("ESTIMATE POWER WITH upstrap, for fixed target effect size = ", eff_tar))
  
  # update sample to represent target effect (here: no update)
  dat_upd <- dat
  dat_upd$y <- dat_upd$y + (eff_tar - fixef(fit_obs)["x1"]) * dat_upd$x1
  # define object to store values across B resamplings
  mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
  ##
  for (bb in 1 : B_boot){ 
    # resample data indices for current boot repetition (upstrap up to N target max)
    dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar_max, replace = TRUE)
    dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar_max, replace = TRUE)
    dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_bb_x1_is0$subjid <- rep(1 : N_tar_max, each = ni) 
    dat_bb_x1_is1$subjid <- rep((N_tar_max + 1) : (N_tar_max + N_tar_max), each = ni)
    dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar_max, each = ni)
    dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar_max, each = ni)
    for (rr in 1 : N_tar_grid_l){ 
      tryCatch({ # bb <- 100; rr <- 1
        N_tar  <- N_tar_grid[rr]
        size_0 <- round(N_tar * 2 * (1-cov_prop))
        size_1 <- (N_tar * 2) - size_0
        dat_bb_x1_is0_bb_rr <- dat_bb_x1_is0 %>% filter(subjid_arm <= size_0)
        dat_bb_x1_is1_bb_rr <- dat_bb_x1_is1 %>% filter(subjid_arm <= size_1)
        # final data set
        dat_bb_rr <- rbind(dat_bb_x1_is0_bb_rr, dat_bb_x1_is1_bb_rr)
        fit_bb_rr  <- lmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr)  
        pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 5]
        mat_boot[rr, bb] <- (pval_bb_rr < 0.05) * 1
      }, error = function(e) {message(e)})
    }
  }
  value <- rowMeans(mat_boot, na.rm = TRUE)
  ##
  mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
  mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
  mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
  mat_out_tmp$name          <- "upstrap_power"
  mat_out_tmp$eff_tru       <- eff_tru
  mat_out_tmp$eff_tar       <- eff_tar
  mat_out_tmp$value         <- value
  mat_out_tmp$sd_sigma      <- sd_sigma
  mat_out_tmp$sd_tau        <- sd_tau
  mat_out_tmp$prop_success  <- mean(!is.na(mat_boot))
  mat_out_tmp$cov_prop      <- cov_prop
  mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
  rm(mat_out_tmp, value, mat_boot)
} 



# ------------------------------------------------------------------------------
# ADDITIONAL SIMULATIONS TO BE RUN ONCE

if (arrayjob_idx == 1){
  message(paste0("RUN ADDITIONAL SIMULATIONS TO BE RUN ONCE"))
  
  for (cov_prop in cov_prop_grid){ # cov_prop  <- cov_prop_grid[2]; bb <- 1
    print(paste0("cov_prop = ", cov_prop))
    
    for (bb in 1 : R_powertrue){
      set.seed(bb)
      x1_i      <- c(rep(0, N_tar_max), rep(1, N_tar_max))
      x2_i      <- rbinom(n = N_tar_max * 2, size = 1, prob = 0.5)
      x3_i      <- runif(n = N_tar_max * 2, min = 0, max = 1)
      b0_i      <- rnorm(n = (N_tar_max * 2), mean = 0, sd = sd_tau)
      subjid_i     <- 1 : (N_tar_max + N_tar_max)      # subject ID unique in data set 
      subjid_arm_i <- c(1 : N_tar_max, 1 : N_tar_max)  # subject ID unique within "treatment arm" 
      # observation-specific
      x1_ij     <- rep(x1_i, each = ni) 
      x2_ij     <- rep(x2_i, each = ni) 
      x3_ij     <- rep(x3_i, each = ni) 
      b0_ij     <- rep(b0_i, each = ni) 
      eps_ij    <- rnorm(n = ((N_tar_max + N_tar_max) * ni), mean = 0, sd = sd_sigma)
      subjid_ij     <- rep(subjid_i, each = ni)            
      subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
      rm(x1_i, x2_i, x3_i, b0_i, subjid_i, subjid_arm_i)
      y_ij      <- (b0_ij + coef_x0) + (eff_tar * x1_ij) + (coef_x2 * x2_ij) +  (coef_x3 * x3_ij) + eps_ij
      dat_bb    <- data.frame(y = y_ij, x1 = x1_ij, x2 = x2_ij,  x3 = x3_ij,
                              subjid = subjid_ij, subjid_arm = subjid_arm_ij)
      value <- rep(NA, N_tar_grid_l)
      for (rr in 1 : N_tar_grid_l){  # rr <- 2
        tryCatch({
          N_tar  <- N_tar_grid[rr]
          size_0 <- round(N_tar * 2 * (1-cov_prop))
          size_1 <- (N_tar * 2) - size_0
          dat_bb_rr_is0 <- dat_bb %>% filter(x1 == 0, subjid_arm <= size_0)
          dat_bb_rr_is1 <- dat_bb %>% filter(x1 == 1, subjid_arm <= size_1)
          dat_bb_rr <- rbind(dat_bb_rr_is0, dat_bb_rr_is1)
          # data set 
          fit_bb_rr  <- lmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr)  
          pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 5]
          value[rr] <- (pval_bb_rr < 0.05) * 1
        }, error = function(e) {message(e)})
      }
      
      mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
      mat_out_tmp$N_obs         <- N_tar_grid
      mat_out_tmp$arrayjob_idx  <- rep(bb, N_tar_grid_l)
      mat_out_tmp$name          <- "true_power"
      mat_out_tmp$eff_tru       <- eff_tar
      mat_out_tmp$eff_tar       <- eff_tar
      mat_out_tmp$value         <- value
      mat_out_tmp$sd_sigma      <- sd_sigma
      mat_out_tmp$sd_tau        <- sd_tau
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
