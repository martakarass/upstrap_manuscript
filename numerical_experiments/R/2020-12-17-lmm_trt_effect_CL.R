#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# print("Started installation of needed packages...")
# install.packages("lme4", repos = "http://cran.us.r-project.org")
# install.packages("lmerTest", repos = "http://cran.us.r-project.org")
# install.packages("here", repos = "http://cran.us.r-project.org")
# print("Installed needed packages.")
print("Started loading needed packages...")
library(lme4)
library(lmerTest)
library(here)
library(matrixStats)
print("Loaded needed packages.")

# read the arguments
arg_str <- as.character(args[1]) # arg_str <- "1"
rep_idx <- as.numeric(arg_str)

# define results dir
out_fpath <- paste0(here::here(), "/numerical_experiments/results_CL/2020-12-17-lmm_trt_effect/repidx_", rep_idx, ".csv")

# define experiment params 
## data 
x_coef  <- 0.5
tau2    <- 1
sigma2  <- 1
n0      <- 50   # sample size of each of the arms
ni      <- 3
## grid
n1_min  <- 50
n1_max  <- 100
# n1_grid   <- n1_min:n1_max
n1_grid   <- seq(from = n1_min, to = n1_max, by = 5)
n1_grid_l <- length(n1_grid)
## upstrap
B_boot <- 100
n0 <- 50  


# ------------------------------------------------------------------------------
# simulate the data for this array job 
set.seed(rep_idx)

subjid_i  <- 1:(n0 * 2)
x_i       <- c(rep(1, n0), rep(0, n0))
x_ij      <- rep(x_i, each = ni) 
subjid_ij <- rep(subjid_i, each = ni) 
b_i       <- rnorm(n = (n0 * 2), mean = 0, sd = tau2)
b_ij      <- rep(b_i, each = ni) 
eps_ij    <- rnorm(n = (n0 * 2 * ni), mean = 0, sd = sigma2)
y_ij      <- b_ij + x_coef * x_ij + eps_ij
dat       <- data.frame(x = x_ij, y = y_ij, subjid = subjid_ij)

# matrix of IDs corrresponding to treatment/non-treatment arm 
dat_subjid_x1 <- 1:50
dat_subjid_x0 <- 51:100

# define object to store simulation results
mat_out <- matrix(NA, nrow = B_boot, ncol = n1_grid_l)

t1 <- Sys.time()
# iterate over boostrap repetitions 
for (B_boot_idx in 1:B_boot){
  if (B_boot_idx %% 10 == 0) print(paste0("B_boot_idx: ", B_boot_idx))
  print(B_boot_idx)
  ## upsample data for current boot repetition
  dat_subjid_x1_b <- sample(x = dat_subjid_x1, size = n1_max, replace = TRUE)
  dat_subjid_x0_b <- sample(x = dat_subjid_x0, size = n1_max, replace = TRUE)
  dat_x1_b <- lapply(dat_subjid_x1_b, function(subjid_tmp) dat[which(subjid_ij %in% subjid_tmp), 1:2]) 
  dat_x0_b <- lapply(dat_subjid_x0_b, function(subjid_tmp) dat[which(subjid_ij %in% subjid_tmp), 1:2]) 
  dat_x1_b <- do.call("rbind", dat_x1_b)
  dat_x0_b <- do.call("rbind", dat_x0_b)
  dat_x1_b$subjid <- rep(1:n1_max, each = ni)
  dat_x0_b$subjid <- rep((n1_max + 1):(2 * n1_max), each = ni)
  dat_x1_b$n1 <- rep(1:n1_max, each = ni)
  dat_x0_b$n1 <- rep(1:n1_max, each = ni)
  dat_b <- rbind(dat_x1_b, dat_x0_b)
  rm(dat_subjid_x1_b, dat_subjid_x0_b, dat_x1_b, dat_x0_b)
  # iterate over n1 values grid
  for (n1_grid_idx in 1:n1_grid_l){ # n1_grid_idx <- 5
    n1_tmp <- n1_grid[n1_grid_idx]
    dat_b_tmp <- dat_b[dat_b$n1 <= n1_tmp, ]
    fit <- lmer(y ~ x + (1 | subjid), data = dat_b_tmp)
    fit_s <- summary(fit)
    x_coef_est <- fit_s$coefficients[2, 5]
    # check if effect identified
    mat_out[B_boot_idx, n1_grid_idx] <- (x_coef_est < 0.05)
  }
  if (B_boot_idx %% 10 == 0){
    # print message 
    secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
    message(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", secs_passed, " secs"))
    # save file (current status)
    power_out <- colMeans(mat_out, na.rm = TRUE)
    power_out <- round(power_out, 5)
    power_df <- data.frame(power_vec = power_out, n1_vec = n1_grid)
    power_df$repidx_vec <- rep_idx
    saveRDS(object = power_df, file = out_fpath)
  }
}

secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
message(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", secs_passed, " secs"))

# save file (current status)
power_out <- colMeans(mat_out, na.rm = TRUE)
power_out <- round(power_out, 5)
power_df <- data.frame(power_vec = power_out, n1_vec = n1_grid)
power_df$repidx_vec <- rep_idx
print(power_df)
saveRDS(object = power_df, file = out_fpath)
print(paste0("SAVED AND FINISHED."))