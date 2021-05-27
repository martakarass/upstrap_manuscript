
rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to data saves
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-twosample_ttest_raw")
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-05-25-twosample_ttest_compare_upstrap_powerttest_agg")
dir.create(res_fdir_agg)

# experiment parameters (relevant subset; see /R/2021-02-17-twosample_ttest.R)
N1_min <- 5
N1_max <- 300 
N1_grid <- N1_min : N1_max

# parameters for bootstrap CI computation
B_bootci     <- 1000
conf_bootci  <- 0.95
alpha_bootci <- (1 + c(-conf_bootci, conf_bootci))/2
# power_val    <- 0.8

# read precomputed raw results data 
fnames_all <- list.files(res_fdir_raw, full.names = TRUE)
fnames_all_1 <- fnames_all[grepl("mat_out_upstrap", fnames_all)]
fnames_all_2 <- sapply(fnames_all_1, function(val) gsub("mat_out_upstrap", "mat_out_powerttest", val))



# ------------------------------------------------------------------------------
# PART 1: bootstrap to get CI for the power estimates: mean, median of the difference

message("Starting PART 1...")
t1_part1 <- Sys.time()

# object to store aggregated results
out_agg_df <- data.frame()

# iterate over data files; aggregate; combine aggregated results
for (i in 1 : length(fnames_all_1)){ # i <- 3
  message(paste0("i: ", i))

  # make data file with differences of estimates
  dat_i_1 <- readRDS(fnames_all_1[i]) # upstrap
  dat_i_2 <- readRDS(fnames_all_2[i]) # power.t.test
  dat_i <- dat_i_1 - dat_i_2
  # summary(as.vector(dat_i))
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  dat_i_method_name <- "diff_ups_powerttest"
  dat_i_namesplit <- strsplit(basename(fnames_all_1[i]), "_")[[1]]
  dat_i_N0 <- as.numeric(gsub(".rds", "", dat_i_namesplit[5]))

  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = dat_i_ncol)
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = dat_i_ncol)
  set.seed(123)
  for (b in 1 : B_bootci){ # b <- i
    if (b %% 100 == 0) message(b)
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
    method_name = rep(dat_i_method_name, dat_i_ncol),
    N0 = rep(dat_i_N0, dat_i_ncol),
    N1 = N1_grid,
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
saveRDS(out_agg_df, paste0(res_fdir_agg, "/twosample_ttest_compare_upstrap_powerttest_bootCI.rds"))
rm(out_agg_df)

message(paste0("Time part 1: ", Sys.time() - t1_part1))
