#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# rm(list = ls()); arg_str <- "1" 

library(lme4)
library(lmerTest)
library(here)
library(matrixStats)
library(dplyr)
library(geepack)
print("Loaded needed packages.")

# read the command line arguments
arg_str <- as.character(args[1]) 
arrayjob_idx <- as.numeric(arg_str)

# define results dir
out_fpath <- paste0(here::here(), "/numerical_experiments/results_CL/2020-12-29-lmm_trt/arrayjob_", arrayjob_idx)

# ------------------------------------------------------------------------------
# define experiment params 
coef_x1 <- 0.5
tau2    <- 1
sigma2  <- 1
N       <- 50   # sample size of each of the two arms
ni      <- 3
N1_min  <- 50
N1_max  <- 100
B_boot  <- 1000 # TODO
innerloop_N <- 100 # TODO

# make object to store simulation results
mat_out <- data.frame()

# set seed for reproducibility
set.seed(arrayjob_idx)

for (innerloop_idx in 1:innerloop_N){
  message(paste0("innerloop_idx: ", innerloop_idx))
  
  # ----------------------------------------------------------------------------
  # simulate the data for this array job, this inner loop 
  
  # below is done for max value of N1 considered to allow to estimate "gold standard" 
  # independent sampling for max N1 considered
  # simulated/generated quantities
  N_tmp <- N1_max
  subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
  subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
  x1_i      <- c(rep(1, N_tmp), rep(0, N_tmp))
  b0_i      <- rnorm(n = (N_tmp * 2), mean = 0, sd = tau2)
  eps_ij    <- rnorm(n = (N_tmp * 2 * ni), mean = 0, sd = sigma2)
  # simulated/generated quantities to data frame vectors 
  subjid_ij     <- rep(subjid_i, each = ni) 
  subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
  x1_ij     <- rep(x1_i, each = ni) 
  b0_ij     <- rep(b0_i, each = ni) 
  y_ij      <- b0_ij + coef_x1 * x1_ij + eps_ij
  dat       <- data.frame(y = y_ij, x1 = x1_ij, subjid = subjid_ij, subjid_arm = subjid_arm_ij)
  
  # define N1 grid
  N1_grid   <- seq(from = N1_min, to = N1_max, by = 3)
  N1_grid_l <- length(N1_grid)
  
  # prepare objects to store results from (a) gold standard, (b) upstrap
  mat_out_tmp               <- data.frame(N0 = rep(N, N1_grid_l), N1 = N1_grid) 
  mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
  mat_out_tmp$innerloop_idx <- rep(innerloop_idx, N1_grid_l)
  
  # ----------------------------------------------------------------------------
  # estimate "gold standard" power 
  
  # define object to store results
  mat_out_gold_LMM   <- matrix(NA, nrow = N1_grid_l, ncol = 1)
  mat_out_gold_GEE   <- matrix(NA, nrow = N1_grid_l, ncol = 1)
  # iterate over number of independent units (subjects) in the same arm
  for (N1_grid_idx in 1:N1_grid_l){ # N1_grid_idx <- 10
    N1_tmp <- N1_grid[N1_grid_idx]
    dat_sub <- dat %>% filter(subjid_arm <= N1_tmp)
    # LMM
    fit_LMM   <- lmer(y ~ x1 + (1 | subjid), data = dat_sub)
    fit_LMM_s <- summary(fit_LMM)
    fit_LMM_coef_est  <- fit_LMM_s$coefficients[2, 5]
    mat_out_gold_LMM[N1_grid_idx, 1] <- (fit_LMM_coef_est < 0.05)
    # GEE
    fit_GEE_formula <- formula(y ~ x1)
    fit_GEE <- geeglm(formula = fit_GEE_formula, family = gaussian(link = "identity"), 
                      data = dat_sub, id = dat_sub$subjid, corstr = "exchangeable")
    fit_GEE_s <- summary(fit_GEE)
    fit_GEE_coef_est  <- fit_GEE_s$coefficients[2, 4]
    mat_out_gold_GEE[N1_grid_idx, 1] <- (fit_GEE_coef_est < 0.05)
  }
  # add results to mat_out
  mat_out_tmp$out_gold_LMM <- mat_out_gold_LMM[, 1] * 1
  mat_out_tmp$out_gold_GEE <- mat_out_gold_GEE[, 1] * 1
  rm(mat_out_gold_LMM, mat_out_gold_GEE)
    
  # ----------------------------------------------------------------------------
  # estimate "upstrap" power 
  
  # data with "observed" sample of size N (base for upstraping) 
  datN <-  dat %>% filter(subjid_arm <= N); rm(dat)
  # matrix of IDs corresponding to treatment/non-treatment arm 
  dat_subjid_x1_is1 <- unique(datN %>% filter(x1 == 1) %>% pull(subjid))
  dat_subjid_x1_is0 <- unique(datN %>% filter(x1 == 0) %>% pull(subjid))
  
  # define object to store results
  mat_out_boot_LMM   <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
  mat_out_boot_GEE   <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
  
  t1 <- Sys.time()
  # iterate over bootstrap repetitions 
  for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 2
    if (B_boot_idx %% 100 == 0){
      secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
      message(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", secs_passed, " secs"))
    }
    ## upsample data for current boot repetition (upstrap up to N1 max)
    dat_subjid_x1_is1_b <- sample(x = dat_subjid_x1_is1, size = N1_max, replace = TRUE)
    dat_subjid_x1_is0_b <- sample(x = dat_subjid_x1_is0, size = N1_max, replace = TRUE)
    dat_x1_is1_b <- lapply(dat_subjid_x1_is1_b, function(subjid_tmp) datN[which(datN$subjid %in% subjid_tmp), 1:2]) 
    dat_x1_is0_b <- lapply(dat_subjid_x1_is0_b, function(subjid_tmp) datN[which(datN$subjid %in% subjid_tmp), 1:2]) 
    dat_x1_is1_b <- do.call("rbind", dat_x1_is1_b)
    dat_x1_is0_b <- do.call("rbind", dat_x1_is0_b)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_x1_is1_b$subjid <- rep(1:N1_max, each = ni)
    dat_x1_is0_b$subjid <- rep((N1_max + 1):(2 * N1_max), each = ni)
    dat_x1_is1_b$subjid_arm <- rep(1:N1_max, each = ni)
    dat_x1_is0_b$subjid_arm <- rep(1:N1_max, each = ni)
    dat_b <- rbind(dat_x1_is1_b, dat_x1_is0_b)
    rm(dat_subjid_x1_is1_b, dat_subjid_x1_is0_b, dat_x1_is1_b, dat_x1_is0_b)
    # iterate over N1 values grid
    for (N1_grid_idx in 1:N1_grid_l){ # N1_grid_idx <- 1
      N1_tmp <- N1_grid[N1_grid_idx]
      dat_sub <- dat_b %>% filter(subjid_arm <= N1_tmp)
      # LMM
      fit_LMM   <- lmer(y ~ x1 + (1 | subjid), data = dat_sub)
      fit_LMM_s <- summary(fit_LMM)
      fit_LMM_coef_est  <- fit_LMM_s$coefficients[2, 5]
      mat_out_boot_LMM[N1_grid_idx, B_boot_idx] <- (fit_LMM_coef_est < 0.05)
      # GEE
      fit_GEE_formula <- formula(y ~ x1)
      fit_GEE <- geeglm(formula = fit_GEE_formula, family = gaussian(link = "identity"), 
                        data = dat_sub, id = dat_sub$subjid, corstr = "exchangeable")
      fit_GEE_s <- summary(fit_GEE)
      fit_GEE_coef_est  <- fit_GEE_s$coefficients[2, 4]
      mat_out_boot_GEE[N1_grid_idx, B_boot_idx] <- (fit_GEE_coef_est < 0.05)
    }
  }
  # add results to mat_out
  mat_out_tmp$out_boot_LMM <- rowMeans(mat_out_boot_LMM, na.rm = TRUE)
  mat_out_tmp$out_boot_GEE <- rowMeans(mat_out_boot_GEE, na.rm = TRUE)
  
  # ----------------------------------------------------------------------------
  # store resutls to master file
  mat_out <- rbind(mat_out, mat_out_tmp)
  rm(mat_out_tmp)
}

# save final result to file 
saveRDS(object = mat_out, file = out_fpath)
print(paste0("SAVED AND FINISHED."))