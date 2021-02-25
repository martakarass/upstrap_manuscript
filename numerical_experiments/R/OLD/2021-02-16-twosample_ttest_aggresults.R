
#' This script aggregates the estimates of power of rejecting H0 in two-sample 
#' t-test problem. It computes aggregates across R repetitions of the experiment:
#' - median + bootstrap CI
#' - mean + bootstrap CI
#' for each sample size considered in the experiment.

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to data saves
res_fdir_agg <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-15-twosample_ttest_agg")
res_fdir_raw <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-15-twosample_ttest_raw")

# experiment parameters (relevant subset)
N1_max <- 250 

# parameters for bootstrap CI computation
B_bootci <- 1000
conf_bootci <- 0.95
alpha_bootci <- (1 + c(-conf_bootci, conf_bootci))/2

# read precomputed raw results data 
fnames_all <- list.files(res_fdir_raw, full.names = TRUE)
dat_all_list <- lapply(fnames_all, readRDS)

# object to store aggregated results
out_agg_df <- data.frame()

# iterate over data files; aggregate; combine aggregated results
for (i in 1 : length(dat_all_list)){ # i <- 1
  message(paste0("i: ", i))
  
  # read i-th file with data 
  dat_i <- dat_all_list[[i]]
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  dat_i_namesplit <- strsplit(basename(fnames_all[i]), "_")[[1]]
  dat_i_method_name <- dat_i_namesplit[3]
  dat_i_N0 <- as.numeric(gsub(".rds", "", dat_i_namesplit[5]))
  dat_i_N1_grid <- dat_i_N0 : N1_max
  
  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = dat_i_nrow)
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = dat_i_nrow)
  set.seed(123)
  for (b in 1 : B_bootci){
    print(b)
    resample_b_idx <- sample(1 : dat_i_nrow, replace = TRUE)
    dat_i_resample_b <- dat_i[resample_b_idx, ]
    dat_i_bootobj_median[b, ] <- matrixStats::colMedians(dat_i_resample_b)
    dat_i_bootobj_mean[b, ]   <- matrixStats::colMeans2(dat_i_resample_b)
  }
  
  # bootstrap CI for median
  dat_i_bootci_median_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[1])
  dat_i_bootci_median_upr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[2])
  dat_i_median <-  matrixStats::colMedians(dat_i)
  # bootstrap CI for mean
  dat_i_bootci_mean_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[1])
  dat_i_bootci_mean_upr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[2])
  dat_i_mean <-  matrixStats::colMeans2(dat_i)
  t2 <- Sys.time()
  
  # put aggregated results into one data frame
  dat_i_agg_df <- data.frame(
    method_name = rep(dat_i_method_name, dat_i_nrow),
    N0 = rep(dat_i_N0, dat_i_nrow),
    N1 = dat_i_N1_grid,
    power_est_aggmedian = dat_i_median,
    power_est_aggmedian_lwr = dat_i_bootci_median_lwr,
    power_est_aggmedian_upr = dat_i_bootci_median_upr,
    power_est_aggmean = dat_i_mean,
    power_est_aggmean_lwr = dat_i_bootci_mean_lwr,
    power_est_aggmean_upr = dat_i_bootci_mean_upr
  )
  # append i-th file-specific aggregated results 
  out_agg_df <- rbind(out_agg_df, dat_i_agg_df)
  
}


# save the aggregated data 
saveRDS(out_agg_df, paste0(res_fdir_agg, "/twosample_ttest_out_agg_df.rds"))
