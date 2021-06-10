#!/usr/bin/env Rscript

#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave 2021-06-08-glmm_trt.R -t 1-1000 -tc 60 -N JOB_glmm_trt
#' qalter -tc 40 job_IDX 

arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)
# rm(list = ls()); arrayjob_idx <- 1

library(lme4)
library(lmerTest)
library(here)
library(matrixStats)
library(dplyr)
library(geepack)
message("Loaded needed packages.")

# define results dir; create dir if any does not exist
out_dir_raw <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-08-glmm_trt_raw")
dir.create(path = out_dir_raw)
message(paste0("dir.exists(path = out_dir_raw): ", dir.exists(path = out_dir_raw)))

# define results path
out_fpath_raw <- paste0(out_dir_raw, "/arrayjob_", arrayjob_idx, ".rds")
message(out_fpath_raw)


# ------------------------------------------------------------------------------
# define experiment params 
coef_0  <- -0.25
coef_x1 <- 0.5
tau2    <- 1
sigma2  <- 1
N0      <- 80   # sample size of each of the two arms
ni      <- 3
N1_min  <- 10
N1_max  <- 250

# define N1 grid
N1_grid   <- seq(from = N1_min, to = N1_max, by = 1)
N1_grid_l <- length(N1_grid)
# number of resamplings 
B_boot    <- 1000 # TODO


# ----------------------------------------------------------------------------
# simulate the data 

# seed seed for reproducibility of the results 
set.seed(arrayjob_idx)

# deterministic quantities
N_tmp <- N1_max
subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
x1_i         <- c(rep(1, N_tmp), rep(0, N_tmp))
# simulated quantities
b0_i      <- rnorm(n = (N_tmp * 2), mean = 0, sd = tau2)
# simulated/generated data frame variables 
subjid_ij     <- rep(subjid_i, each = ni) 
subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
x1_ij     <- rep(x1_i, each = ni) 
b0_ij     <- rep(b0_i, each = ni) 
XB_ij     <- coef_0 + b0_ij + coef_x1 * x1_ij
p_ij      <- 1/(1 + exp(-XB_ij))
y_ij      <- rbinom(n = length(p_ij), size = 1, prob = p_ij)
dat       <- data.frame(y = y_ij, x1 = x1_ij, subjid = subjid_ij, subjid_arm = subjid_arm_ij)

# make object to store simulation results 
mat_out <- data.frame()

# make object to keep track of the beginning of time 
t1 <- Sys.time()



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATION: UPSTRAP POWER

message(paste0("--- Running: uptstrap N0=", N0))

# data with "observed" sample of size N0 
datN0 <-  dat %>% filter(subjid_arm <= N0)
# subj_id corresponding to treatment/non-treatment arm 
dat_subjid_x1_is1 <- unique(datN0 %>% filter(x1 == 1) %>% pull(subjid))
dat_subjid_x1_is0 <- unique(datN0 %>% filter(x1 == 0) %>% pull(subjid))

# define object to store results
# mat_out_boot_LMM   <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
mat_out_boot_GEE   <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)

# iterate over bootstrap repetitions 
for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 2
  print(paste0("B_boot_idx: ", B_boot_idx))
  if (B_boot_idx %% 10 == 0){
    t_passed <- round(as.numeric(Sys.time() - t1, unit = "mins"))
    message(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", t_passed, " mins"))
  }
  ## upsample data for current boot repetition (upstrap up to N1 max)
  dat_subjid_x1_is1_b <- sample(x = dat_subjid_x1_is1, size = N1_max, replace = TRUE)
  dat_subjid_x1_is0_b <- sample(x = dat_subjid_x1_is0, size = N1_max, replace = TRUE)
  dat_x1_is1_b <- lapply(dat_subjid_x1_is1_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:2]) 
  dat_x1_is0_b <- lapply(dat_subjid_x1_is0_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:2]) 
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
  for (N1_grid_idx in 1:N1_grid_l){ # N1_grid_idx <- 200
    tryCatch({
      dat_b_sub <- dat_b[dat_b$subjid_arm <= N1_grid[N1_grid_idx], ]
      # GEE
      fit_GEE_formula <- formula(y ~ x1)
      fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                        data = dat_b_sub, id = dat_b_sub$subjid, corstr = "exchangeable")
      fit_GEE_s <- summary(fit_GEE)
      fit_GEE_coef_est  <- fit_GEE_s$coefficients[2, 4]
      mat_out_boot_GEE[N1_grid_idx, B_boot_idx] <- (fit_GEE_coef_est < 0.05) * 1
      rm(dat_b_sub)
    }, error = function(e) {message(e)})
  }
  rm(dat_b)
}

# store and append the results 
mat_out_tmp               <- data.frame(N0 = rep(N0, N1_grid_l), N1 = N1_grid) 
mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
mat_out_tmp$name          <- "upstrap_power_GEE"
mat_out_tmp$value         <- rowMeans(mat_out_boot_GEE, na.rm = TRUE)

mat_out <- rbind(mat_out, mat_out_tmp); rm(mat_out_tmp, mat_out_boot_GEE)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATION: BOOTSTRAP POWER

message(paste0("--- Running: bootstrap"))

# define object to store values across B resamplings
mat_out_boot_GEE <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)

for (N1_grid_idx in 1 : length(N1_grid)){ # N1_grid_idx <-1
  if (N1_grid_idx %% 10 == 0){
    message(paste0("N1_grid_idx: ", N1_grid_idx))
  }
  N1 <- N1_grid[N1_grid_idx]
  N0 <- N1  # bootstrap
  # data with "observed" sample of size N0 
  datN0 <-  dat %>% filter(subjid_arm <= N0)
  # subj_id corresponding to treatment/non-treatment arm 
  dat_subjid_x1_is1 <- unique(datN0 %>% filter(x1 == 1) %>% pull(subjid))
  dat_subjid_x1_is0 <- unique(datN0 %>% filter(x1 == 0) %>% pull(subjid))
  
  # iterate over bootstrap repetitions 
  for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 10
    ## upsample data for current boot repetition (upstrap up to N1 max)
    dat_subjid_x1_is1_b <- sample(x = dat_subjid_x1_is1, size = N1, replace = TRUE)
    dat_subjid_x1_is0_b <- sample(x = dat_subjid_x1_is0, size = N1, replace = TRUE)
    dat_x1_is1_b <- lapply(dat_subjid_x1_is1_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:2]) 
    dat_x1_is0_b <- lapply(dat_subjid_x1_is0_b, function(subjid_tmp) datN0[which(datN0$subjid %in% subjid_tmp), 1:2]) 
    dat_x1_is1_b <- do.call("rbind", dat_x1_is1_b)
    dat_x1_is0_b <- do.call("rbind", dat_x1_is0_b)
    ## make new subj ID so as to treat resampled subjects as new ones
    dat_x1_is1_b$subjid <- rep(1:N1, each = ni)
    dat_x1_is0_b$subjid <- rep((N1 + 1):(2 * N1), each = ni)
    dat_x1_is1_b$subjid_arm <- rep(1:N1, each = ni)
    dat_x1_is0_b$subjid_arm <- rep(1:N1, each = ni)
    dat_b <- rbind(dat_x1_is1_b, dat_x1_is0_b)
    tryCatch({
      # # LMM
      # fit_LMM   <- lmer(y ~ x1 + (1 | subjid), data = dat_b)
      # fit_LMM_s <- summary(fit_LMM)
      # fit_LMM_coef_est  <- fit_LMM_s$coefficients[2, 5]
      # mat_out_boot_LMM[N1_grid_idx, B_boot_idx] <- (fit_LMM_coef_est < 0.05)
      # GEE
      fit_GEE_formula <- formula(y ~ x1)
      fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                        data = dat_b, id = dat_b$subjid, corstr = "exchangeable")
      fit_GEE_s <- summary(fit_GEE)
      fit_GEE_coef_est  <- fit_GEE_s$coefficients[2, 4]
      mat_out_boot_GEE[N1_grid_idx, B_boot_idx] <- (fit_GEE_coef_est < 0.05) * 1
    }, error = function(e) {message(e)})
    rm(dat_b)
  }
}

# store and append the results 
mat_out_tmp               <- data.frame(N0 = rep(N0, N1_grid_l), N1 = N1_grid) 
mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
mat_out_tmp$name          <- "bootstrap_power_GEE"
mat_out_tmp$value         <- rowMeans(mat_out_boot_GEE, na.rm = TRUE)

mat_out <- rbind(mat_out, mat_out_tmp); rm(mat_out_tmp, mat_out_boot_GEE)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATION: SINGLE RESULT 

message(paste0("--- Running: test result"))

# define object to store values across B resamplings
mat_out_boot_GEE <- matrix(NA, nrow = N1_grid_l, ncol = 1)

for (N1_grid_idx in 1 : length(N1_grid)){ # N1_grid_idx <-1
  if (N1_grid_idx %% 10 == 0){
    message(paste0("N1_grid_idx: ", N1_grid_idx))
  }
  N1 <- N1_grid[N1_grid_idx]
  N0 <- N1 
  # data with "observed" sample of size N0 
  datN0 <-  dat %>% filter(subjid_arm <= N0)
  # run one result
  fit_GEE_formula <- formula(y ~ x1)
  fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                    data = datN0, id = datN0$subjid, corstr = "exchangeable")
  fit_GEE_s <- summary(fit_GEE)
  fit_GEE_coef_est  <- fit_GEE_s$coefficients[2, 4]
  mat_out_boot_GEE[N1_grid_idx, 1] <- (fit_GEE_coef_est < 0.05) * 1
  rm(datN0)
}

# add results to mat_out
mat_out_tmp               <- data.frame(N0 = N1_grid, N1 = N1_grid) 
mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
mat_out_tmp$name          <- "run_result"
mat_out_tmp$value         <- rowMeans(mat_out_boot_GEE, na.rm = TRUE)

# store results to master file
mat_out <- rbind(mat_out, mat_out_tmp); rm(mat_out_tmp)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SAVE FINAL RESULT TO FILE
saveRDS(object = mat_out, file = out_fpath_raw)
t_diff <- round(difftime(Sys.time(), t1, units = "mins"), 3)
message(paste0("Saved FINAL state. Minutes elapsed: ", t_diff))

# Saved FINAL state. Minutes elapsed: 3.296


