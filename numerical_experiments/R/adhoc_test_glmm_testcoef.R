#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' rm $ups/numerical_experiments/results_CL_shared/2021-08-04-glmm_adhoc_test/*
#' cd $ups/numerical_experiments/R
#' ls -l -d *adhoc*
#' rm JOB_adhoc_glmm*
#' 
#' Rnosave adhoc_test_glmm_testcoef.R -t 1-50 -tc 50 -N JOB_adhoc_glmm


arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# rm(list = ls()); arrayjob_idx <- 1
set.seed(arrayjob_idx)

library(lme4)
library(lmerTest)
library(here)
library(tidyverse)
library(simr)

dir_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-04-glmm_adhoc_test")
dir.create(dir_out)

# experiment parameters
N_obs   <- 50 # number of subjects in each of the two "treatment arms" 
ni      <- 3  # number of observations per subject 
coef_x0 <- -0.5 
coef_x1 <- 0.6
coef_x2 <- 0.1 
coef_x3 <- -0.1
tau2    <- 1
sigma2  <- 1

N_tar   <- 50
eff_tar <- 0.6

# number of boot repetitions within one experiment, one setup
B_boot  <- 500
R_rep   <- 20

result_glmm      <- rep(NA, R_rep)
power_upstrap_v1 <- rep(NA, R_rep)
power_upstrap_v2 <- rep(NA, R_rep)
power_upstrap_v3 <- rep(NA, R_rep)
power_upstrap_v4 <- rep(NA, R_rep)
power_simr       <- rep(NA, R_rep)


# ------------------------------------------------------------------------------

t1 <- Sys.time()

for (rr_rep in 1 : R_rep){
  
  message(paste0("REPETITION: ", rr_rep, " -- ",  round(as.numeric(Sys.time() - t1, unit = "mins"), 1), " mins"))
  
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
  # SINGLE RUN
  result_glmm[rr_rep] <- as.numeric(summary(fit_obs)$coef["x1", 4] < 0.05)
  
  
  # ------------------------------------------------------------------------------
  # [ v1 ]
  # ESTIMATE POWER WITH upstrap, for fixed target effect size 
  
  message(paste0("upstrap -- v1"))
  
  dat_upd <- dat
  dat_upd$link_orig <- predict(fit_obs, type = "link") 
  dat_upd$link_upd  <- dat_upd$link_orig  + (eff_tar - fixef(fit_obs)["x1"]) * dat_upd$x1
  dat_upd$res_upd   <- 1/(1 + exp(-dat_upd$link_upd))
  mat_boot <- numeric( B_boot)
  for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
    if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    # resample data indices for current boot repetition (upstrap up to N target max)
    dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar, replace = TRUE)
    dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar, replace = TRUE)
    dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat_upd %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_bb_x1_is1$subjid <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid <- rep((N_tar + 1) : (2 * N_tar), each = ni)
    dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar, each = ni)
    # length(unique(c(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjid)))
    # intersect(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjidy)
    dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
    # simulate response on resampled data, assuming target effect size  
    dat_bb$y <- rbinom(n = nrow(dat_bb), size = 1, prob = dat_bb$res_upd)
    # iterate over target sample (here: only one)
    dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
    fit_bb_rr  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
    pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
    mat_boot[bb] <- (pval_bb_rr < 0.05) * 1
  }
  power_upstrap_v1[rr_rep] <- mean(mat_boot)
  
  
  # ------------------------------------------------------------------------------
  # [ v2 ]
  # ESTIMATE POWER WITH upstrap, for fixed target effect size 
  
  # message(paste0("upstrap -- v2"))
  # 
  # rm(fit_obs)
  # mat_boot <- numeric( B_boot)
  # for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
  #   if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
  #   # resample data indices for current boot repetition (upstrap up to N target max)
  #   dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar, replace = TRUE)
  #   dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar, replace = TRUE)
  #   dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
  #   dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
  #   ## make new subj ID so as to treat resampled subjects as new ones
  #   dat_bb_x1_is1$subjid <- rep(1 : N_tar, each = ni)
  #   dat_bb_x1_is0$subjid <- rep((N_tar + 1) : (2 * N_tar), each = ni)
  #   dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar, each = ni)
  #   dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar, each = ni)
  #   # length(unique(c(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjid)))
  #   # intersect(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjidy)
  #   dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
  #   fit_bb <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb, family = binomial)  
  #   dat_bb$link_orig <- predict(fit_bb, type = "link") 
  #   dat_bb$link_upd  <- dat_bb$link_orig  + (eff_tar - fixef(fit_bb)["x1"]) * dat_bb$x1
  #   dat_bb$res_upd   <- 1/(1 + exp(-dat_bb$link_upd))
  #   dat_bb$y         <- rbinom(n = nrow(dat_bb), size = 1, prob = dat_bb$res_upd)
  #   # iterate over target sample (here: only one)
  #   dat_bb_rr    <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
  #   fit_bb_rr    <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
  #   pval_bb_rr   <- summary(fit_bb_rr)$coef["x1", 4]
  #   mat_boot[bb] <- (pval_bb_rr < 0.05) * 1
  #   # summary(fit_bb)
  #   # summary(fit_bb_rr)
  # }
  # power_upstrap_v2[rr_rep] <- mean(mat_boot)
  
  
  # ------------------------------------------------------------------------------
  # [ v3 ]
  # ESTIMATE POWER WITH upstrap, for fixed target effect size 
  
  message(paste0("upstrap -- v3"))

  rm(fit_obs)
  mat_boot <- numeric( B_boot)
  for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
    if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    # resample data indices for current boot repetition (upstrap up to N target max)
    dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar, replace = TRUE)
    dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar, replace = TRUE)
    dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_bb_x1_is1$subjid <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid <- rep((N_tar + 1) : (2 * N_tar), each = ni)
    dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar, each = ni)
    dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
    # length(unique(c(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjid)))
    # intersect(dat_bb_x1_is1$subjid, dat_bb_x1_is0$subjidy)
    #ee
    fit_bb <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb, family = binomial)  
    fixef(fit_bb)["x1"] <- eff_tar
    # update y
    dat_bb$y   <- simulate(fit_bb, nsim = 1)[[1]]
    # subset target sample (starts from max)
    dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
    fit_bb_rr  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
    pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
    mat_boot[bb] <- (pval_bb_rr < 0.05) * 1
  }
  power_upstrap_v3[rr_rep] <- mean(mat_boot)
  
  
  # ------------------------------------------------------------------------------
  # [ v-4 ]
  # ESTIMATE POWER WITH upstrap, for fixed target effect size 
  
  fit_obs <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat, family = binomial)  
  fit_obs_upd <- fit_obs
  fixef(fit_obs_upd)["x1"] <- eff_tar
  mat_boot <- numeric( B_boot)
  for (bb in 1 : B_boot){ # bb <- 1; rr <- 1
    if (bb %% 100 == 0) message(paste0("bb: ", bb, " [", round(bb / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    # resample data indices for current boot repetition (upstrap up to N target max)
    dat_bb_idx_x1_is1 <- sample(x = dat_subjid_x1_is1, size = N_tar, replace = TRUE)
    dat_bb_idx_x1_is0 <- sample(x = dat_subjid_x1_is0, size = N_tar, replace = TRUE)
    dat_bb_x1_is1     <- lapply(dat_bb_idx_x1_is1, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    dat_bb_x1_is0     <- lapply(dat_bb_idx_x1_is0, function(subjid_tmp) dat %>% filter(subjid == subjid_tmp)) %>% do.call("rbind", .)
    # simulate new Y 
    dat_bb_x1_is1$y <- simulate(fit_obs_upd, nsim = 1, newdata = dat_bb_x1_is1)[[1]]
    dat_bb_x1_is0$y <- simulate(fit_obs_upd, nsim = 1, newdata = dat_bb_x1_is0)[[1]]
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_bb_x1_is1$subjid <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid <- rep((N_tar + 1) : (2 * N_tar), each = ni)
    dat_bb_x1_is1$subjid_arm <- rep(1 : N_tar, each = ni)
    dat_bb_x1_is0$subjid_arm <- rep(1 : N_tar, each = ni)
    dat_bb <- rbind(dat_bb_x1_is1, dat_bb_x1_is0)
    # subset target sample (starts from max)
    dat_bb_rr  <- dat_bb[dat_bb$subjid_arm <= N_tar, ]
    fit_bb_rr  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat_bb_rr, family = binomial)  
    pval_bb_rr <- summary(fit_bb_rr)$coef["x1", 4]
    mat_boot[bb] <- (pval_bb_rr < 0.05) * 1
  }
  power_upstrap_v4[rr_rep] <- mean(mat_boot)
  
  # ------------------------------------------------------------------------------
  # ESTIMATE POWER WITH simr, for fixed target effect size 
  
  message(paste0("SIMR"))
  
  fit_obs  <- glmer(y ~ x1 + x2 + x3 + (1 | subjid), data = dat, family = binomial)
  fit_obs_simr <- fit_obs
  # update effect size
  fixef(fit_obs_simr)["x1"] <- eff_tar
  ps_out <- simr::powerSim(fit_obs_simr, nsim = B_boot, progress = FALSE)
  ps_out_s      <- summary(ps_out)
  power_simr[rr_rep] <- ps_out_s$mean
  
  
  # ------------------------------------------------------------------------------
  # SAVE TO FILE 
  df <- data.frame(result_glmm = result_glmm, 
                   power_upstrap_v1 = power_upstrap_v1, 
                   power_upstrap_v2 = power_upstrap_v2,
                   power_upstrap_v3 = power_upstrap_v3,
                   power_upstrap_v4 = power_upstrap_v4,
                   power_simr = power_simr)
  df$arrayjob_idx <- arrayjob_idx
  df$lanuch <- 3
  saveRDS(object = df, file = paste0(dir_out, "/arrayjob_", arrayjob_idx, ".rds"))
  
  rm(fit_obs)
}

