#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave run_glmm_testcoef.R -t 1-1000 -tc 100 -N JOB_glmm


arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# rm(list = ls()); arrayjob_idx <- 1
set.seed(arrayjob_idx)
message(paste0("arrayjob_idx: ", arrayjob_idx))

library(lme4)
library(lmerTest)
library(here)
library(tidyverse)
library(matrixStats)
library(simr)

# dir to save results 
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-02-glmm_testcoef_raw")
# create dirs if any does not exist
dir.create(path = res_fdir_raw)

# experiment parameters
N_obs   <- 50 # number of subjects in each of the two "treatment arms" 
ni      <- 3  # number of observations per subject 
coef_x0 <- -0.5 
coef_x1 <- 0.6
coef_x2 <- 0.1 
coef_x3 <- -0.1
tau2    <- 1
sigma2  <- 1
N_tar_min     <- 10
N_tar_max     <- 150
N_tar_grid    <- seq(N_tar_min, N_tar_max, by = 3)
N_tar_grid_l  <- length(N_tar_grid)
eff_tar_grid  <- c(0.4, 0.6, 0.8)
eff_tru       <- coef_x1

# number of boot repetitions within one experiment, one setup
B_boot  <- 1000

# simulate sample (for maximum sample size first)
# subject-specific
x1_i      <- rep(c(0, 1), times = N_obs)
x2_i      <- rbinom(n = N_obs * 2, size = 1, prob = 0.5)
x3_i      <- runif(n = N_obs * 2, min = 0, max = 1)
b0_i      <- rnorm(n = (N_obs * 2), mean = 0, sd = tau2)
subjid_i     <- 1 : (N_obs * 2)           # subject ID unique in data set 
subjid_arm_i <- rep(1 : N_obs, each = 2)  # subject ID unique within "treatment arm" 
# observation-specific
x1_ij     <- rep(x1_i, each = ni) 
x2_ij     <- rep(x2_i, each = ni) 
x3_ij     <- rep(x3_i, each = ni) 
b0_ij     <- rep(b0_i, each = ni) 
subjid_ij     <- rep(subjid_i, each = ni)            
subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
XB_ij     <- (b0_ij + coef_x0) + (coef_x1 * x1_ij) + (coef_x2 * x2_ij) +  (coef_x3 * x3_ij)
p_ij      <- 1/(1 + exp(-XB_ij))
y_ij      <- rbinom(n = length(p_ij), size = 1, prob = p_ij)
# data set 
dat       <- data.frame(y = y_ij, x1 = x1_ij, x2 = x2_ij,  x3 = x3_ij,
                        subjid = subjid_ij, subjid_arm = subjid_arm_ij)

# indexes of individuals in each two "treatment arms" 
dat_subjid_x1_is1 <- unique(dat %>% filter(x1 == 1) %>% pull(subjid))
dat_subjid_x1_is0 <- unique(dat %>% filter(x1 == 0) %>% pull(subjid))

# get observed effect 
fit_obs  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat, family = binomial)

message(paste0("nrow(dat): ", nrow(dat)))
message(paste0("mean(dat$y): ", round(mean(dat$y), 2)))
message(paste0("fixef(fit_obs)['x1']: ", round(fixef(fit_obs)["x1"], 2)))


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
for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
  print(bb)
  if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
  # resample data indices for current boot repetition (upstrap up to N target max)
  dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar_max, replace = TRUE)
  dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar_max, replace = TRUE)
  dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
  dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
  ## make new subj ID so as to treat resampled subjects as new ones
  dat_bb_x1_is1$subjid <- rep(1 : N_tar_max, each = ni)
  dat_bb_x1_is0$subjid <- rep((N_tar_max + 1) : (2 * N_tar_max), each = ni)
  dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar_max, each = ni)
  dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar_max, each = ni)
  # length(unique(c(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjid)))
  # intersect(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjidy)
  dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
  for (rr in 1 : N_tar_grid_l){ 
    tryCatch({
      N_tar      <- N_tar_grid[rr]
      dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
      fit_bb_rr  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
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
  # define link term value
  dat_upd$link_orig <- predict(fit_obs, type = "link") 
  # define link term value for  target effect size
  dat_upd$link_upd  <- dat_upd$link_orig  + (eff_tar - fixef(fit_obs)["x1"]) * dat_upd$x1
  dat_upd$res_upd   <- 1/(1 + exp(-dat_upd$link_upd))
  # define object to store values across B resamplings
  mat_boot <- matrix(NA, nrow = N_tar_grid_l, ncol = B_boot)
  ##
  for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
    print(bb)
    if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    # resample data indices for current boot repetition (upstrap up to N target max)
    dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar_max, replace = TRUE)
    dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar_max, replace = TRUE)
    dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_bb_x1_is1$subjid <- rep(1 : N_tar_max, each = ni)
    dat_bb_x1_is0$subjid <- rep((N_tar_max + 1) : (2 * N_tar_max), each = ni)
    dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar_max, each = ni)
    dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar_max, each = ni)
    # length(unique(c(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjid)))
    # intersect(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjidy)
    dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
    # simulate response on resampled data, assuming target effect size  
    dat_bb$y <- rbinom(n = nrow(dat_bb), size = 1, prob = dat_bb$res_upd)
    for (rr in 1 : N_tar_grid_l){ 
      tryCatch({
        N_tar      <- N_tar_grid[rr]
        dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
        fit_bb_rr  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
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
  mat_out_tmp$eff_tar       <- eff_tar # observed power case
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
  fixef(fit_obs_simr)["x1"] <- eff_tar
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


# table(mat_out_all$name, mat_out_all$eff_tar, useNA = 'always')

t2 <- Sys.time()
t_diff <- round(as.numeric(t2 - t1, unit = "mins"), 1)
message(paste0("COMPLETED. Time difference of ", t_diff, " minutes."))
# Time difference of 8 minutes.


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SAVE TO FILE 

out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
saveRDS(object = mat_out_all, file = out_fpath_raw)