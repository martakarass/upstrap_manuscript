
rm(list = ls())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' @description 
#' Script to compare sample size needed in two 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>


## -----------------------------------------------------------------------------

## args
# rm(list = ls())
# arg_str <- as.character(args[1])
arg_str <- "N0_200_MU_01_R_50_BBOOT_10000"

## fixed params
# project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
project_dir <- "/users/mkaras/_PROJECTS/upstrap_manuscript"

# derivative args
n0 <- as.numeric(strsplit(arg_str, split = "_")[[1]][2])
R  <- as.numeric(strsplit(arg_str, split = "_")[[1]][6])
B_boot  <- as.numeric(strsplit(arg_str, split = "_")[[1]][8])
# function to inset characters at certain index of a string 
insert_chars_at_index <- function(old, index, chars) {
  paste(substr(old, 1, index-1), as.character(chars), substr(old, index, nchar(old)), sep = "")
}
mu0 <- strsplit(arg_str, split = "_")[[1]][4]
mu  <- as.numeric(insert_chars_at_index(mu0, 2, "."))

# file name to save the results to 
out_fname <- paste0(arg_str, ".csv")
out_fpath <- paste0(project_dir,"/numerical_experiments/results_CL/2020-10-09-ttest_normal_mean_results/", out_fname)
# message with derivative args values
message(paste0(c("N0", "mu", "R", "B_boot: "), collapse = ", "))
message(paste0(c(n0, mu, R, B_boot), collapse = ", "))

# libraries
library(data.table)
library(dplyr)



## -----------------------------------------------------------------------------
## FUNCTIONS

# function to get approximation of correction factor c4
# https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
get_c4_n <- function(n){
  c4_n <- 1 - (1/(4 * n)) -  (7/(32 * n^2)) - (19/(128 * n^3))
  return(c4_n)
}

# get bias-corrected estimator of sample standard deviation based on sd
get_ubiased_sd_est <- function(sd_s, n){
  c4_n <- get_c4_n(n)
  out <- sd_s/c4_n
  return(out)
}

# function to compute t.test (~ 3x faster than t.test() function)
t.test.pval_FAST <- function(x){
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  df <- nx-1
  stderr <- sqrt(vx/nx)
  tstat <- mx/stderr
  pval <- 2 * pt(-abs(tstat), df)
  return(pval)
}


## -----------------------------------------------------------------------------

# power (fixed param)
power_val <- 0.8
set.seed(1)

rr <- 1
sample_obs <- rnorm(n = n0, mean = mu, sd = 1)
delta_obs  <- mean(sample_obs)
sd_obs     <- sd(sample_obs)

## Approach 1: power.t.test
n1_powerttest_out <- power.t.test(delta = delta_obs, sd = sd_obs, power = power_val)
n1_powerttest <- ceiling(n1_powerttest_out$n)

## Approach 2: upstrap
n1_min <- 100
n1_max <- n1_powerttest * 10
boot_resamples <- matrix(sample(x = sample_obs, size = (B_boot * n1_max), replace = TRUE), nrow = B_boot, ncol = n1_max)

find_n1_RECURSIVE <- function(resamples_tmp, n1_min_tmp, n1_max_tmp, tol_lvl = 0.001){
  n1_tmp <- round(mean(c(n1_min_tmp, n1_max_tmp)))
  boot_resamples_tmp_pval  <- apply(resamples_tmp[, 1:n1_tmp], 1, t.test.pval_FAST)
  boot_resamples_tmp_power <- mean(boot_resamples_tmp_pval < 0.05)
  print(paste0("current params: n in [", n1_min_tmp, ",", n1_max_tmp, "] tmp: ", n1_tmp, " power: ", round(boot_resamples_tmp_power, 6)))
  print(abs(boot_resamples_tmp_power - power_val))
  if ((abs(boot_resamples_tmp_power - power_val) < tol_lvl) | ((n1_max_tmp - n1_min_tmp) < 2)) {
    return(n1_tmp)
  } else if (boot_resamples_tmp_power < power_val){
    return(find_n1_RECURSIVE(resamples_tmp, n1_tmp, n1_max_tmp, tol_lvl))
  } else {
    return(find_n1_RECURSIVE(resamples_tmp, n1_min_tmp, n1_tmp, tol_lvl))
  }
}
  
find_n1_RECURSIVE(boot_resamples, n1_min, n1_max)





## ---------------------------------------------------------------------------

set.seed(2)
boot_resamples_TRUE <- matrix(rnorm(10000 * 10000, mean = mu), nrow = 10000, ncol = 10000)
find_n1_RECURSIVE(boot_resamples_TRUE, n1_min, n1_max, tol_lvl = 0.0001)






