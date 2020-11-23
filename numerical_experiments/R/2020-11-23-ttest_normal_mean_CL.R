#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' @description 
#' Script to bootstrap-estimate standard deviation / variance of estimator "s" 
#' of standard deviation. 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>


## -----------------------------------------------------------------------------


# project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
project_dir <- "/users/mkaras/_PROJECTS/upstrap_manuscript"
results_dir <- paste0(project_dir, "/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean")
library(data.table)

# params 
arg_str <- as.character(args[1])
# arg_str <- "MU_0.1_N0_200_R_1000_BBOOT_10000_N1MIN_200_N1MAX_1200_REPIDX1_99001"
# var vals provided
mu        <- as.numeric(strsplit(arg_str, split = "_")[[1]][2])
n0        <- as.numeric(strsplit(arg_str, split = "_")[[1]][4])
rep_n     <- as.numeric(strsplit(arg_str, split = "_")[[1]][6])
B_boot    <- as.numeric(strsplit(arg_str, split = "_")[[1]][8])
n1_min    <- as.numeric(strsplit(arg_str, split = "_")[[1]][10])
n1_max    <- as.numeric(strsplit(arg_str, split = "_")[[1]][12])
rep_idx_1 <- as.numeric(strsplit(arg_str, split = "_")[[1]][14])
# var vals derivative
n1_grid     <- n1_min:n1_max
n1_grid_max <- max(n1_grid)
# print var vals 
print(paste0(
  "mu: ", mu, 
  ", n0: ", n0, 
  ", rep_n: ", rep_n, 
  ", B_boot: ", B_boot, 
  ", n1_min: ", n1_min, 
  ", n1_max: ", n1_max,
  ", rep_idx_1: ", rep_idx_1))

# function to compute cumulative var
cumvar <- function (x, sd = FALSE) {
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  # if (sd) v <- sqrt(v)
  v
}

# compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
vals_cum_reject_H0 <- function(vals){
  vals_cumnx      <- seq_along(vals)
  vals_cumdf      <- vals_cumnx - 1
  vals_cummean    <- cumsum(vals) / vals_cumnx
  vals_cumvar     <- cumvar(vals)
  vals_cumstderr  <- sqrt(vals_cumvar / vals_cumnx)
  vals_cumststat  <- vals_cummean / vals_cumstderr
  vals_cumpval    <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
  vals_reject_H0  <- vals_cumpval < 0.05
  return(vals_reject_H0)
}

# define objects to store simulation results
mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = length(n1_grid))
mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = length(n1_grid))
print("Dimmensions of objects defined to store simulation results")
dim(mat_out_powerttest)
dim(mat_out_upstrap)

# define paths to save simulation results
mat_out_powerttest_fpath <- paste0(results_dir, "/out_powerttest_repidx1_", rep_idx_1, ".csv")
mat_out_upstrap_fpath    <- paste0(results_dir, "/out_upstrap_repidx1_", rep_idx_1, ".csv")

print("Starting the simulation...")
t1 <- Sys.time()
for (i in 1:rep_n){
  # set seed unique across all jobs 
  rep_idx_i <- rep_idx_1 + (i - 1)
  set.seed(rep_idx_i)
  # simulate sample 
  sample_i <- rnorm(n = n0, mean = mu, sd = 1)
  # sample_i: generate power estimate with power.t.test()
  sample_i_mean <- mean(sample_i)
  sample_i_sd   <- sd(sample_i)
  mat_out_powerttest[i, ] <- sapply(n1_grid, function(n_tmp){
    (power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
  })
  mat_out_powerttest[i, ] <- round(mat_out_powerttest[i, ], 5)
  # sample_i: generate power estimate with upstrap()
  # generate the samples
  boot_resamples_i <- matrix(sample(x = sample_i, size = (B_boot * n1_grid_max), replace = TRUE), 
                             nrow = B_boot, ncol = n1_grid_max)
  boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
  mat_out_upstrap[i, ] <- apply(boot_resamples_i_rejectH0[, n1_grid], 2, mean)
  mat_out_upstrap[i, ] <- round(mat_out_upstrap[i, ], 5)
  # save to file every 100
  if (i %% 50 == 0){
    fwrite(as.data.table(mat_out_powerttest), mat_out_powerttest_fpath)
    fwrite(as.data.table(mat_out_upstrap),    mat_out_upstrap_fpath)
    secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
    paste0("i: ", i, " [", round(i / rep_n * 100, 2), "%], ", secs_passed, " secs")
  }
  
}
t2 <- Sys.time()

# final save
fwrite(as.data.table(mat_out_powerttest), mat_out_powerttest_fpath)
fwrite(as.data.table(mat_out_upstrap),    mat_out_upstrap_fpath)
secs_passed <- round(as.numeric(t2 - t1, unit = "secs"))
paste0("i: ", i, " [", round(i / rep_n * 100, 2), "%], ", secs_passed, " secs")
print("SIMULATION COMPLETED.")





