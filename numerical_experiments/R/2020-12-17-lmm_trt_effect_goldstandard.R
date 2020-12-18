
rm(list = ls())

library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(matrixStats)
out_fpath <- paste0(here::here(), "/numerical_experiments/results/2020-12-17-lmm_trt_effect/power_df_goldstandard.rda")

x_coef <- 0.5
tau2   <- 1
sigma2 <- 1
n0     <- 50   # sample size of each of the arms
ni     <- 3

rep_n  <- 100 * 100
n1_min    <- 50
n1_max    <- 100
n1_grid   <- n1_min:n1_max
n1_grid_l <- length(n1_grid)

mat_out_goldstandard <- matrix(NA, nrow = rep_n, ncol = length(n1_grid))

set.seed(1)
t1 <- Sys.time()
for (n1_tmp_idx in 1:n1_grid_l){ 
  n1_tmp    <- n1_grid[n1_tmp_idx]
  subjid_i  <- 1:(n1_tmp * 2)
  x_i       <- c(rep(1, n1_tmp), rep(0, n1_tmp))
  x_ij      <- rep(x_i, each = ni) 
  subjid_ij <- rep(subjid_i, each = ni) 
  for (i in 1:rep_n){
    if (i %% 100 == 0) print(paste0("i: ", i))
    b_i    <- rnorm(n = (n1_tmp * 2), mean = 0, sd = tau2)
    b_ij   <- rep(b_i, each = ni) 
    eps_ij <- rnorm(n = (n1_tmp * 2 * ni), mean = 0, sd = sigma2)
    y_ij   <-  b_ij + x_coef * x_ij + eps_ij
    dat    <-  data.frame(x = x_ij, y = y_ij, subjid = subjid_ij)
    # fit the model
    fit <- lmer(y_ij ~ x_ij + (1 | subjid_ij))
    fit_s <- summary(fit)
    x_coef_est <- fit_s$coefficients[2, 5]
    # check if effect identified
    mat_out_goldstandard[i, n1_tmp_idx] <- (x_coef_est < 0.05)
  }
  # print message
  secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
  message(paste0("n1_tmp_idx: ", n1_tmp_idx, " [", round(n1_tmp_idx / n1_grid_l * 100, 2), "%], ", secs_passed, " secs"))
}
t2 <- Sys.time()
round(as.numeric(t2 - t1, unit = "secs")) # 252

power_out_goldstandard <- colMeans(mat_out_goldstandard, na.rm = TRUE)
power_out_goldstandard <- round(power_out_goldstandard, 5)
power_df <- data.frame(power_vec = power_out_goldstandard, n1_vec = n1_grid)
saveRDS(object = power_df, file = out_fpath)

