
#' This script aggregates the estimates of power of rejecting H0 in LMM
#' regression problem.
#' 
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave 2021-02-22-lmm_trt_aggresults.R -N JB_lmm_trt_agg


rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to data saves
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-lmm_trt_raw")
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-lmm_trt_agg")
dir.create(path = res_fdir_agg)

# experiment parameters (relevant subset; see /R/2021-02-17-lmm_trt.R)
N  <- 41   # sample size of each of the two arms

# parameters for bootstrap CI computation
B_bootci     <- 1000
conf_bootci  <- 0.95
alpha_bootci <- (1 + c(-conf_bootci, conf_bootci))/2
power_val    <- 0.8

# read precomputed raw results data 
fnames_all   <- list.files(res_fdir_raw, full.names = TRUE)
dat_all <- do.call("rbind", lapply(fnames_all, readRDS))
dat_all <- mutate(dat_all, idx = paste0(arrayjob_idx, "_", innerloop_idx))     
dim(dat_all)
head(dat_all)
str(dat_all)
names(dat_all)

# how many out of R=1000 has processed
dat_all %>% select(arrayjob_idx, innerloop_idx)  %>% distinct() %>% nrow()
length(unique(dat_all$idx))
# Mar 1: 960


# ------------------------------------------------------------------------------
# PART 1: bootstrap to get CI for the power estimates: mean, median

message("Starting PART 1...")

ylab_vec <- c("out_gold_LMM", "out_gold_GEE", "out_boot_LMM", "out_boot_GEE")
out_df_all <- data.frame()

t1 <- Sys.time()
for (i in 1 : length(ylab_vec)){ # ylab_i <- "out_boot_LMM"
  message(paste0("PART 1: outer i: ", i))
  # current outcome
  ylab_i <- ylab_vec[i]

  # matrix with current outcome [number of repetitions R] x [number of unique values of N1]
  dat_i <- 
    dat_all %>% 
    select(all_of(c('idx', 'N1', ylab_i))) %>% 
    rename_at(vars(ylab_i), function(x) "y") %>% 
    pivot_wider(names_from = N1, values_from = y) %>% 
    select(-idx) %>%
    as.matrix()
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  
  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = dat_i_ncol)
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = dat_i_ncol)
  set.seed(123)
  for (b in 1 : B_bootci){
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
  # prepare final data frame with these results
  out_df_i <- data.frame(
    est_median = dat_i_median,
    est_bootci_median_lwr = dat_i_bootci_median_lwr,
    est_bootci_median_upr = dat_i_bootci_median_upr,
    est_mean = dat_i_mean,
    est_bootci_mean_lwr = dat_i_bootci_mean_lwr,
    est_bootci_mean_upr = dat_i_bootci_mean_upr
  ) %>% mutate(
    method_name = rep(ylab_i, dat_i_ncol),
    N0 = rep(N, dat_i_ncol),
    N1 = as.numeric(colnames(dat_i)), 
    .before = everything(),
  )
  # append results
  out_df_all <- rbind(out_df_all, out_df_i)
}
t2 <- Sys.time()
t2 - t1
# Time difference of 51.94879 secs

# save the aggregated data 
saveRDS(out_df_all, paste0(res_fdir_agg, "/lmm_trt_power_est_bootCI.rds"))
rm(out_df_all)


# ------------------------------------------------------------------------------
# PART 2: bootstrap minimum sample size needed to get certain power met 

message("Starting PART 2...")

ylab_vec <- c("out_gold_LMM", "out_gold_GEE", "out_boot_LMM", "out_boot_GEE")
out_df_all <- data.frame()
idx_unique <- unique(dat_all$idx)

t1 <- Sys.time()
for (i in 1 : length(ylab_vec)){ # ylab_i <- "out_boot_LMM"
  message(paste0("PART 2: outer i: ", i))
  # current outcome
  ylab_i <- ylab_vec[i]

  # matrix with current outcome [number of repetitions R] x [number of unique values of N1]
  dat_i <- 
    dat_all %>% 
    select(all_of(c('idx', 'N1',  ylab_i))) %>% 
    rename_at(vars(ylab_i), function(x) "y")  %>%
    pivot_wider(names_from = N1, values_from = y) %>% 
    select(-idx) %>%
    as.matrix()
  dat_i_nrow <- nrow(dat_i)
  dat_i_ncol <- ncol(dat_i)
  
  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = 1)
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = 1)
  set.seed(123)
  for (b in 1 : B_bootci){ # b <- 1
    if (b %% 100 == 0) message(b)
    resample_b_idx <- sample(1 : dat_i_nrow, replace = TRUE)
    # resample from wide, make a new unique row index `idx`, reshape to long
    dat_i_resample_b <- 
      dat_i[resample_b_idx, ] %>%
      as.data.frame() %>%
      mutate(idx = row_number()) %>%
      pivot_longer(cols = -idx, names_to = "N1", values_to = "y") %>%
      mutate(N1 = as.numeric(N1))
    dat_i_bootobj_median[b, 1] <- dat_i_resample_b %>% group_by(N1) %>% summarize(y = median(y)) %>% filter(y >= power_val) %>% pull(N1) %>% min()
    dat_i_bootobj_mean[b, 1]   <- dat_i_resample_b %>% group_by(N1) %>% summarize(y = mean(y))   %>% filter(y >= power_val) %>% pull(N1) %>% min()
    # dat_i_bootobj_median[b, 1] <- dat_i_resample_b %>% group_by(idx) %>% filter(y >= power_val) %>% filter(N1 == min(N1)) %>% pull(N1) %>% median()
    # dat_i_bootobj_mean[b, 1]   <- dat_i_resample_b %>% group_by(idx) %>% filter(y >= power_val) %>% filter(N1 == min(N1)) %>% pull(N1) %>% mean()
  }
  # bootstrap CI for median
  dat_i_bootci_median_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[1])
  dat_i_bootci_median_upr <-  matrixStats::colQuantiles(dat_i_bootobj_median, probs = alpha_bootci[2])
  # bootstrap CI for mean
  dat_i_bootci_mean_lwr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[1])
  dat_i_bootci_mean_upr <-  matrixStats::colQuantiles(dat_i_bootobj_mean, probs = alpha_bootci[2])
  # estimates: median, mean
  dat_i_long <- 
    dat_all %>% 
    select(all_of(c('idx', 'N1',  ylab_i))) %>% 
    rename_at(vars(ylab_i), function(x) "y")
  dat_i_median <- dat_i_long %>% group_by(N1) %>% summarize(y = median(y)) %>% filter(y >= power_val) %>% pull(N1) %>% min()
  dat_i_mean   <- dat_i_long %>% group_by(N1) %>% summarize(y = mean(y))   %>% filter(y >= power_val) %>% pull(N1) %>% min()
  
  # prepare final data frame with these results
  out_df_i <- data.frame(
    est_median = dat_i_median,
    est_bootci_median_lwr = dat_i_bootci_median_lwr,
    est_bootci_median_upr = dat_i_bootci_median_upr,
    est_mean = dat_i_mean,
    est_bootci_mean_lwr = dat_i_bootci_mean_lwr,
    est_bootci_mean_upr = dat_i_bootci_mean_upr
  ) %>% mutate(
    method_name = rep(ylab_i, 1),
    N0 = rep(N, 1),
    N1 = NA, 
    .before = everything(),
  )
  # append results
  out_df_all <- rbind(out_df_all, out_df_i)
}
t2 <- Sys.time()
t2 - t1
# Time difference of 51.94879 secs

# save the aggregated data 
saveRDS(out_df_all, paste0(res_fdir_agg, "/lmm_trt_samplesize_est_bootCI.rds"))
rm(out_df_all)
