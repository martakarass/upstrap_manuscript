
#' This script aggregates the estimates of power of rejecting H0 in one-sample
#' t-test problem.
#' 
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave 2021-02-22-onesample_ttest_aggresults.R -N JB_onesample_agg

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to data saves
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-onesample_ttest_raw")
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-onesample_ttest_agg")

# experiment parameters (relevant subset; see /R/2021-02-17-onesample_ttest.R)
N1_min <- 5
N1_max <- 200 
N1_grid <- N1_min : N1_max

# parameters for bootstrap CI computation
B_bootci     <- 1000
conf_bootci  <- 0.95
alpha_bootci <- (1 + c(-conf_bootci, conf_bootci))/2
power_val    <- 0.8

# read precomputed raw results data 
fnames_all <- list.files(res_fdir_raw, full.names = TRUE)
dat_all_list <- lapply(fnames_all, readRDS)


# ------------------------------------------------------------------------------
# PART 1: bootstrap to get CI for the power estimates: mean, median

message("Starting PART 1...")
t1_part1 <- Sys.time()

# object to store aggregated results
out_agg_df <- data.frame()

# iterate over data files; aggregate; combine aggregated results
for (i in 1 : length(dat_all_list)){ # i <- 3
  message(paste0("i: ", i))

  # read i-th file with data 
  dat_i <- dat_all_list[[i]]
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  dat_i_namesplit <- strsplit(basename(fnames_all[i]), "_")[[1]]
  dat_i_method_name <- dat_i_namesplit[3]
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
saveRDS(out_agg_df, paste0(res_fdir_agg, "/onesample_ttest_power_est_bootCI.rds"))
rm(out_agg_df)

message(paste0("Time part 1: ", Sys.time() - t1_part1))

# ------------------------------------------------------------------------------
# PART 2: bootstrap minimum sample size needed to get certain power met 

message("Starting PART 2...")
t1_part2 <- Sys.time()

# object to store aggregated results
out_agg_df <- data.frame()

# iterate over data files; aggregate; combine aggregated results
for (i in 1 : length(dat_all_list)){ # i <- 3
  message(paste0("i: ", i))
  
  # read i-th file with data 
  dat_i <- dat_all_list[[i]]
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  dat_i_namesplit <- strsplit(basename(fnames_all[i]), "_")[[1]]
  dat_i_method_name <- dat_i_namesplit[3]
  dat_i_N0 <- as.numeric(gsub(".rds", "", dat_i_namesplit[5]))
  
  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = 1)
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = 1)
  set.seed(123)
  for (b in 1 : B_bootci){ # b <- i
    if (b %% 100 == 0) message(b)
    resample_b_idx <- sample(1 : dat_i_nrow, replace = TRUE)
    dat_i_resample_b <- 
      dat_i[resample_b_idx, ] %>%
      as.data.frame() %>%
      rename_all(~as.character(N1_grid)) %>%
      mutate(idx = row_number())  %>%
      pivot_longer(cols = -idx, names_to = "N1", values_to = "y") %>%
      mutate(N1 = as.numeric(N1))
    dat_i_bootobj_median[b, 1] <- dat_i_resample_b %>% group_by(N1) %>% summarize(y = median(y)) %>% filter(y >= power_val) %>% pull(N1) %>% min()
    dat_i_bootobj_mean[b, 1]   <- dat_i_resample_b %>% group_by(N1) %>% summarize(y = mean(y))   %>% filter(y >= power_val) %>% pull(N1) %>% min()
  }
  # bootstrap CI for median
  dat_i_bootci_median_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[1])
  dat_i_bootci_median_upr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[2])
  # bootstrap CI for mean
  dat_i_bootci_mean_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[1])
  dat_i_bootci_mean_upr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[2])
  # estimates: median, mean
  dat_i_long <- 
    dat_i  %>%
    as.data.frame() %>%
    rename_all(~as.character(N1_grid)) %>%
    mutate(idx = row_number())  %>%
    pivot_longer(cols = -idx, names_to = "N1", values_to = "y") %>%
    mutate(N1 = as.numeric(N1))
  dat_i_median <- dat_i_long %>% group_by(N1) %>% summarize(y = median(y)) %>% filter(y >= power_val) %>% pull(N1) %>% min()
  dat_i_mean   <- dat_i_long %>% group_by(N1) %>% summarize(y = mean(y))   %>% filter(y >= power_val) %>% pull(N1) %>% min()
  
  # put aggregated results into one data frame
  dat_i_agg_df <- data.frame(
    method_name = rep(dat_i_method_name, 1),
    N0 = rep(dat_i_N0, 1),
    N1 = NA,
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
message(paste0("Time part 2: ", Sys.time() - t1_part2))
message(paste0("Time TOTAL: ", Sys.time() - t1_part1))

saveRDS(out_agg_df, paste0(res_fdir_agg, "/onesample_ttest_samplesize_est_bootCI.rds"))
rm(out_agg_df)

