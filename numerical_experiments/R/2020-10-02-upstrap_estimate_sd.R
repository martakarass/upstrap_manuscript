
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' @description 
#' Script to bootstrap-estimate standard deviation / variance of estimator "s" 
#' of standard deviation. 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>


## -----------------------------------------------------------------------------

## args
# rm(list = ls())
arg_str <- as.character(args[1])
# arg_str <- "N0_100_N1_150_R_50_BBOOT_10000"
# arg_str <- "N0_1000_N1_10000_R_100000_BBOOT_10000"

## fixed params
# project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
project_dir <- "/users/mkaras/_PROJECTS/upstrap_manuscript"

# derivative args
n0 <- strsplit(arg_str, split = "_")[[1]][2]
n1 <- strsplit(arg_str, split = "_")[[1]][4]
R  <- strsplit(arg_str, split = "_")[[1]][6]
B_boot  <- strsplit(arg_str, split = "_")[[1]][8]

n0 <- as.numeric(n0)
n1 <- as.numeric(n1)
R  <- as.numeric(R)
B_boot  <- as.numeric(B_boot)
out_fname <- paste0(arg_str, ".csv")
out_fpath <- paste0(project_dir,"/numerical_experiments/results_CL/2020-10-02-upstrap_estimate_sd_results/", out_fname)
# message with derivative args values
message(paste0(c("N0", "N1", "R", "B_boot"), collapse = ", "))
message(paste0(c(n0, n1, R, B_boot), collapse = ", "))

# libraries
library(data.table)
library(dplyr)


## -----------------------------------------------------------------------------

# function to get approximation of correction factor c4
# https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
get_c4_n <- function(n){
  c4_n <- 1 - (1/(4 * n)) -  (7/(32 * n^2)) - (19/(128 * n^3))
  return(c4_n)
}
get_sd_s <- function(n, sigma = 1){
  c4_n <- get_c4_n(n)
  sd_s <- sigma * sqrt(1 - c4_n^2)
  return(sd_s)
}
# function to compute E of "s" 
get_mean_s <- function(n, sigma = 1){
  c4_n <- get_c4_n(n)
  mean_s <- sigma * c4_n
  return(mean_s)
}

# function to generate ready-to-save summary of current state of simulation
get_out <- function(out_R_boot_stat_means, out_R_boot_stat_sds, rr, t_start){
  
  # generate time summary
  t_tmp <- Sys.time()
  t_diff <- round(as.numeric(difftime(t_tmp, t_start, units = "secs")))
  
  # quantiles' probs for summary of results from R runs 
  probs_vec   <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  # labels for quantiles' probs for summary of results from R runs 
  probs_label <- gsub(pattern = "\\.", replacement = "", as.character(probs_vec, length.out = 5))
  
  # get the summary of results from R runs: out_R_boot_stat_mean
  out_R_boot_stat_means_summary <- c(
    quantile(out_R_boot_stat_means, probs = probs_vec, na.rm = TRUE),
    mean(out_R_boot_stat_means, na.rm = TRUE),
    sd(out_R_boot_stat_means, na.rm = TRUE))
  # get the summary of results from R runs 
  out_R_boot_stat_sds_summary <- c(
    quantile(out_R_boot_stat_sds, probs = probs_vec, na.rm = TRUE),
    mean(out_R_boot_stat_sds, na.rm = TRUE),
    sd(out_R_boot_stat_sds, na.rm = TRUE))
  
  # define final data frame
  out_mat <- c(
    n0, n1, R, B_boot,
    t_diff, rr,
    out_R_boot_stat_means_summary, 
    out_R_boot_stat_sds_summary)
  out_df <- data.frame(t(matrix(out_mat)))
  names(out_df) <- c(
    'n0', 'n1', 'r', 'B_boot',
    't_diff', 'rr',
    paste0('s_mean_rrepsummary_', probs_label),
    paste0('s_mean_rrepsummary_', c("mean", "sd")),
    paste0('s_sd_rrepsummary_', probs_label),
    paste0('s_sd_rrepsummary_', c("mean", "sd"))
  )
  
  return(out_df)
}

# objects to store simulation results
out_R_boot_stat_means <- rep(NA, R)
out_R_boot_stat_sds   <- rep(NA, R)

# perform R repetitions of the experiment
set.seed(1)
t_start <- Sys.time()
for (rr in 1:R){ # rr <- 1
  # simulate observed sample
  vals <- rnorm(n = n0, mean = 0, sd = 1)
  # generate "B_boot" number of resamples
  boot_resamples <- matrix(sample(x = vals, size = (B_boot * n1), replace = TRUE), nrow = B_boot, ncol = n1)
  # for each boot resamples, compute resample statistic (here: sample standard deviation, "s")
  boot_resamples_stat <- apply(boot_resamples, 1, sd)
  # compute mean and sd of "B_boot" number of resample statistics
  out_R_boot_stat_means[rr] <- mean(boot_resamples_stat)
  out_R_boot_stat_sds[rr]   <- sd(boot_resamples_stat)
  # save partial results - at every 1000-th repetition
  if (rr %% 1000 == 0){
    message(paste0("rr: ", rr))
    out_df <- get_out (out_R_boot_stat_means, out_R_boot_stat_sds, rr, t_start)
    fwrite(as.data.table(out_df), out_fpath)
  }
}

# save completed results at the end 
out_df <- get_out (out_R_boot_stat_means, out_R_boot_stat_sds, rr, t_start)
fwrite(as.data.table(out_df), out_fpath)

