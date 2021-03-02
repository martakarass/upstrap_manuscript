#!/usr/bin/env Rscript

# rm(list = ls()); arg_str <- "1"

library(here)
library(matrixStats)
library(tidyverse)
message("Loaded needed packages.")

# read the command line arguments
arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# define results dir
out_dir_raw <- paste0(here::here(), "/numerical_experiments/results_CL/2021-03-01-lm_trt_raw")
# create dirs if any does not exist
dir.create(path = out_dir_raw)
message(paste0("dir.exists(path = out_dir_raw): ", dir.exists(path = out_dir_raw)))
out_fpath_raw <- paste0(out_dir_raw, "/arrayjob_", arrayjob_idx, ".rds")
message(out_fpath_raw)


# ------------------------------------------------------------------------------
# define experiment params 
coef_x0 <- 0 # set to zero 
coef_x1 <- 0.5
coef_x2 <- 0 # set to zero 
coef_x3 <- 0 # set to zero 
sigma2  <- 1
N0_grid <- c(20, 31, 49)   # sample size of each of the two arms
N1_min  <- 10
N1_max  <- 150
# 1         0.3     20
# 2         0.5     31
# 3         0.7     49
# 4         0.9     91
# 5         0.95   113

# define N1 grid
N1_grid   <- seq(from = N1_min, to = N1_max, by = 1)
N1_grid_l <- length(N1_grid)

innerloop_N <- 20 # TODO => -t 1-50
B_boot      <- 1000 # TODO

# make object to store simulation results (specific to this array job)
mat_out <- data.frame()

t1 <- Sys.time()
# set seed for reproducibility
set.seed(arrayjob_idx)

for (innerloop_idx in 1 : innerloop_N){ # innerloop_idx <- 1
  for (N0 in N0_grid){
    message(paste0("innerloop_idx: ", innerloop_idx, ", N0: ", N0))
  
    # ----------------------------------------------------------------------------
    # make object to store simulation results (specific to this array job, this inner loop)
    mat_out_tmp               <- data.frame(N0 = rep(N0, N1_grid_l), N1 = N1_grid) 
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$innerloop_idx <- rep(innerloop_idx, N1_grid_l)
    mat_out_tmp$out_gold_LM   <- NA
  
    # ----------------------------------------------------------------------------
    # simulate the data (specific to this array job, this inner loop) 
    
    # deterministic quantities
    N_tmp <- N1_max
    subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
    subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
    x1_i         <- c(rep(1, N_tmp), rep(0, N_tmp))
    # simulated quantities
    x2_i         <- rbinom(n = N_tmp * 2, size = 1, prob = 0.5)
    x3_i         <- runif(n = N_tmp * 2, min = 18, max = 100)
    eps_i        <- rnorm(N_tmp * 2, sd = sqrt(sigma2))
    # simulated/generated data frame variables 
    y_i   <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) + (coef_x3 * x3_i) + eps_i
    dat   <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i,
                        subjid = subjid_i, subjid_arm = subjid_arm_i)
  
    # ----------------------------------------------------------------------------
    # estimate "gold standard" power 
    
    # iterate over number of independent units (subjects) in the same arm
    for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 10
      tryCatch({
        N1_tmp <- N1_grid[N1_grid_idx]
        dat_sub <- dat %>% filter(subjid_arm <= N1_tmp)
        # LM
        fit_LM   <- lm(y ~ x1 + x2 + x3, data = dat_sub)
        fit_LM_s <- summary(fit_LM)
        fit_LM_coefpval  <- fit_LM_s$coefficients[2, 4]
        mat_out_tmp$out_gold_LM[N1_grid_idx] <- (fit_LM_coefpval < 0.05) * 1
      }, error = function(e) {message(e)})
    }
  
      
    # ----------------------------------------------------------------------------
    # estimate "upstrap" power 
    
    # data with "observed" sample of size N (base for upstraping) 
    datN0 <-  dat %>% filter(subjid_arm <= N0); rm(dat)
    # matrix of IDs corresponding to treatment/non-treatment arm 
    dat_subjid_x1_is1 <- unique(datN0 %>% filter(x1 == 1) %>% pull(subjid))
    dat_subjid_x1_is0 <- unique(datN0 %>% filter(x1 == 0) %>% pull(subjid))
    
    # define object to store results
    mat_out_boot_LM   <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
  
    # iterate over bootstrap repetitions 
    for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 2
      if (B_boot_idx %% 100 == 0){
        secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
        message(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", secs_passed, " secs"))
      }
      ## upsample data for current boot repetition (upstrap up to N1 max)
      dat_subjid_x1_is1_b <- sample(x = dat_subjid_x1_is1, size = N1_max, replace = TRUE)
      dat_subjid_x1_is0_b <- sample(x = dat_subjid_x1_is0, size = N1_max, replace = TRUE)
      dat_x1_is1_b <- lapply(dat_subjid_x1_is1_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:4]) 
      dat_x1_is0_b <- lapply(dat_subjid_x1_is0_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:4]) 
      dat_x1_is1_b <- do.call("rbind", dat_x1_is1_b)
      dat_x1_is0_b <- do.call("rbind", dat_x1_is0_b)
      ## make new subj ID so as to treat resampled subjects as new ones
      dat_x1_is1_b$subjid <- rep(1 : N1_max, each = 1)
      dat_x1_is0_b$subjid <- rep((N1_max + 1) : (2 * N1_max), each = 1)
      dat_x1_is1_b$subjid_arm <- rep(1 : N1_max, each = 1)
      dat_x1_is0_b$subjid_arm <- rep(1 : N1_max, each = 1)
      dat_b <- rbind(dat_x1_is1_b, dat_x1_is0_b)
      rm(dat_subjid_x1_is1_b, dat_subjid_x1_is0_b, dat_x1_is1_b, dat_x1_is0_b)
      # iterate over N1 values grid
      for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 1
        tryCatch({
          N1_tmp <- N1_grid[N1_grid_idx]
          dat_sub <- dat_b %>% filter(subjid_arm <= N1_tmp)
          # LMM
          # LM
          fit_LM   <- lm(y ~ x1 + x2 + x3, data = dat_sub)
          fit_LM_s <- summary(fit_LM)
          fit_LM_coefpval  <- fit_LM_s$coefficients[2, 4]
          mat_out_boot_LM[N1_grid_idx, B_boot_idx] <- (fit_LM_coefpval < 0.05) * 1
        }, error = function(e) {message(e)})
      }
    }
    # add results to mat_out
    mat_out_tmp$out_boot_LM <- rowMeans(mat_out_boot_LM, na.rm = TRUE)
  
    # ----------------------------------------------------------------------------
    # store results to master file
    mat_out <- rbind(mat_out, mat_out_tmp)
    rm(mat_out_tmp)
    
    # save current results
    saveRDS(object = mat_out, file = out_fpath_raw)
    t_diff <- round(difftime(Sys.time(), t1, units = "mins"), 3)
    message(paste0("Saved current state. Minutes elapsed: ", t_diff))
  }
}

# save final result to file 
saveRDS(object = mat_out, file = out_fpath_raw)
t_diff <- round(difftime(Sys.time(), t1, units = "mins"), 3)
message(paste0("Saved FINAL state. Minutes elapsed: ", t_diff))





