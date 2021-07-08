#!/usr/bin/env Rscript

#' This script estimate power of rejecting H0 in one-sample t-test problem. 
#' 
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave run_lm_testcoef.R -t 1-1000 -tc 50 -N JOB_lm


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
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-lm_testcoef_raw")
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-lm_testcoef_agg")
# remove dirs if exist to make a room for new ones
unlink(res_fdir_agg, recursive = TRUE, force = TRUE)
unlink(res_fdir_raw, recursive = TRUE, force = TRUE)
# create dirs if any does not exist
dir.create(path = res_fdir_agg)
dir.create(path = res_fdir_raw)
message(paste0("dir.exists(path = res_fdir_agg): ", dir.exists(path = res_fdir_agg)))
message(paste0("dir.exists(path = res_fdir_raw): ", dir.exists(path = res_fdir_raw)))

# experiment parameters
N_obs   <- 40 # size of each of the two "treatment arms" 
coef_x0 <- 0 
coef_x1 <- 0.6
coef_x2 <- 0.3 
coef_x3 <- -0.1
sigma2  <- 1
N_tar_min     <- 10
N_tar_max     <- 150
N_tar_grid    <- seq(N_tar_min, N_tar_max, by = 1)
N_tar_grid_l  <- length(N_tar_grid)
eff_tar_grid  <- c(0.3, 0.6, 0.9)
eff_tru       <- coef_x1

# number of boot repetitions within one experiment, one setup
B_boot  <- 1000

# simulate sample (for maximum sample size first)
subjid_arm_i <- rep(1 : N_obs, each = 2)
subjid_i  <- 1 : (N_obs * 2)   # subject ID unique in data set 
x1_i      <- rep(c(0, 1), times = N_obs)
x2_i      <- rbinom(n = N_obs * 2, size = 1, prob = 0.5)
x3_i      <- runif(n = N_obs * 2, min = 18, max = 100)
eps_i     <- rnorm(N_obs * 2, sd = sqrt(sigma2))
y_i       <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) +  (coef_x3 * x3_i) + eps_i
dat       <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i, 
                        subjid = subjid_i, subjid_arm = subjid_arm_i)
dim(dat)

# indexes of individuals in each two "treatment arms" 
dat_0_idx <- which(dat$x1 == 0)
dat_1_idx <- which(dat$x1 == 1)

# get observed effect 
fit_obs  <- lm(y ~ x1 + x2 + x3, data = dat)
coef(fit_obs)["x1"]
# 0.7626672 



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATIONS

t1 <- Sys.time()
mat_out_all <- data.frame()

  
# ------------------------------------------------------------------------------
# ESTIMATE POWER WITH upstrap, for observed effect size 

# update sample to represent target effect (here: no update)
dat_upd <- dat
# define object to store values across B resamplings
mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
##
for (bb in 1 : B_boot){ # bb <- 100; rr <- 1
  # print(bb)
  if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
  # resample data indices for current boot repetition (upstrap up to N target max)
  dat_bb_idx <- c(sample(dat_1_idx, size = N_tar_max, replace = TRUE), 
                  sample(dat_0_idx, size = N_tar_max, replace = TRUE))
  dat_bb            <- dat_upd[dat_bb_idx, ]
  dat_bb$subjid     <- 1 : (2 * N_tar_max)                # subject's ID -- unique in data set
  dat_bb$subjid_arm <- c(1 : N_tar_max, 1 : N_tar_max)    # subject's ID -- unique within a "treatment arm"
  for (rr in 1 : N_tar_grid_l){ 
    tryCatch({
      N_tar      <- N_tar_grid[rr]
      dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
      fit_bb_rr  <- lm(y ~ x1 + x2 + x3, data = dat_bb_rr)
      pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
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
mat_out_tmp$eff_tar       <- NA # observed power case
mat_out_tmp$value         <- value
mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
rm(mat_out_tmp, value, mat_boot)



# ------------------------------------------------------------------------------
# ESTIMATE POWER WITH upstrap, for fixed target effect size 

for (eff_tar in eff_tar_grid){ # eff_tar <- eff_tar_grid[1]
  message(paste0("eff_tar = ", eff_tar))
  # update sample to represent target effect (here: no update)
  dat_upd <- dat
  dat_upd$y <- dat_upd$y + (eff_tar - coef(fit_obs)["x1"]) * dat_upd$x1
  # define object to store values across B resamplings
  mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
  ##
  for (bb in 1:B_boot){ # bb <- 100
    # print(bb)
    if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    # resample data for current boot repetition (upstrap up to N target max)
    dat_bb_idx <- c(sample(dat_1_idx, size = N_tar_max, replace = TRUE), 
                    sample(dat_0_idx, size = N_tar_max, replace = TRUE))
    dat_bb            <- dat_upd[dat_bb_idx, ]
    dat_bb$subjid     <- 1 : (2 * N_tar_max)                # subject's ID -- unique in data set
    dat_bb$subjid_arm <- c(1 : N_tar_max, 1 : N_tar_max)    # subject's ID -- unique within a "treatment arm"
    for (rr in 1 : N_tar_grid_l){ 
      tryCatch({
        N_tar      <- N_tar_grid[rr]
        dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
        fit_bb_rr  <- lm(y ~ x1 + x2 + x3, data = dat_bb_rr)
        pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
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
  mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
  rm(mat_out_tmp, value, mat_boot)
}



# ------------------------------------------------------------------------------
# ESTIMATE POWER WITH simr, for observed effect size 

# define simr-specific grid of sample size 
# (represents total number of subjects, not: number subjects in each of the "treatment arms")
N_tar_grid_simr <- N_tar_grid * 2 

fit_obs_simr  <- fit_obs
fit_ext_simr  <- simr::extend(fit_obs_simr, along = "subjid", n = max(N_tar_grid_simr))
pc_out        <- simr::powerCurve(fit_ext_simr, along = "subjid", breaks = N_tar_grid_simr, nsim = B_boot, progress = FALSE)
pc_out_s      <- summary(pc_out)
value         <- pc_out_s$mean
##
mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
mat_out_tmp$name          <- "simr_power"
mat_out_tmp$eff_tru       <- eff_tru
mat_out_tmp$eff_tar       <- NA # observed power case
mat_out_tmp$value         <- value
mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
rm(mat_out_tmp, value, fit_ext, pc_out, pc_out_s)



# ------------------------------------------------------------------------------
# ESTIMATE POWER WITH simr, for fixed target effect size 

for (eff_tar in eff_tar_grid){ 
  message(paste0("SIMR -- eff_tar = ", eff_tar))
  fit_obs_simr <- fit_obs
  # update effect size
  coef(fit_obs_simr)["x1"] <- eff_tar
  fit_ext_simr  <- simr::extend(fit_obs_simr, along = "subjid", n = max(N_tar_grid_simr))
  pc_out        <- simr::powerCurve(fit_ext_simr, along = "subjid", breaks = N_tar_grid_simr, nsim = B_boot, progress = FALSE)
  pc_out_s      <- summary(pc_out)
  value         <- pc_out_s$mean
  ##
  mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
  mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
  mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
  mat_out_tmp$name          <- "simr_power"
  mat_out_tmp$eff_tru       <- eff_tru
  mat_out_tmp$eff_tar       <- eff_tar
  mat_out_tmp$value         <- value
  mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
  rm(mat_out_tmp, value, fit_obs_simr, fit_ext_simr, pc_out, pc_out_s)
}

# t2 <- Sys.time()
# t2 - t1
# Time difference of 2.720334 mins

# table(mat_out_all$name, mat_out_all$eff_tar, useNA = 'always')



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ADDITIONAL SIMULATIONS TO BE RUN ONCE 

if (arrayjob_idx == 1){

  for (eff_tar in eff_tar_grid){ 
    mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
    message(paste0("INDEPSAMPLES -- eff_tar = ", eff_tar))
    for (bb in 1 : B_boot){ # bb <- 100
      if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
      subjid_arm_i <- rep(1 : N_tar_max, each = 2)
      subjid_i  <- 1 : (N_tar_max * 2)   # subject ID unique in data set
      x1_i      <- rep(c(0, 1), times = N_tar_max)
      x2_i      <- rbinom(n = N_tar_max * 2, size = 1, prob = 0.5)
      x3_i      <- runif(n = N_tar_max * 2, min = 18, max = 100)
      eps_i     <- rnorm(N_tar_max * 2, sd = sqrt(sigma2))
      y_i       <- coef_x0 + (eff_tar * x1_i) + (coef_x2 * x2_i) +  (coef_x3 * x3_i) + eps_i
      dat_N_tar_max <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i, 
                                  subjid = subjid_i, subjid_arm = subjid_arm_i)
      for (rr in 1 : N_tar_grid_l){  # rr <- 10
        tryCatch({
          N_tar      <- N_tar_grid[rr]
          dat_bb_rr  <- dat_N_tar_max[dat_N_tar_max$subjid_arm <= N_tar, ]
          fit_bb_rr  <- lm(y ~ x1 + x2 + x3, data = dat_bb_rr)
          pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
          mat_boot[rr, bb] <- (pval_bb_rr < 0.05) * 1
        }, error = function(e) {message(e)})
      }
    }
    value <- rowMeans(mat_boot, na.rm = TRUE)
    ##
    mat_out_tmp               <- data.frame(N_tar = N_tar_grid)
    mat_out_tmp$N_obs         <- rep(N_obs, N_tar_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N_tar_grid_l)
    mat_out_tmp$name          <- "indepsample_power"
    mat_out_tmp$eff_tru       <- eff_tar
    mat_out_tmp$eff_tar       <- eff_tar
    mat_out_tmp$value         <- value
    mat_out_all               <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value, mat_boot)
  }
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SAVE TO FILE 

out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
saveRDS(object = mat_out_all, file = out_fpath_raw)
