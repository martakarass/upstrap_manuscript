
#' This script aggregates the estimates of power of rejecting H0 in LM
#' regression problem.
#' 
#' Notes: 
#' cd $ups 
#' cd numerical_experiments/R
#' Rnosave 2021-06-10-glmm_trt_aggresults.R -N JOB_glmm_agg


rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to precomputed results
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-08-glmm_trt_raw")

# dir to data saves
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-06-08-glmm_trt_agg")
dir.create(path = res_fdir_agg)

# parameters for bootstrap CI computation
B_bootci     <- 1000  # TODO
conf_bootci  <- 0.95
alpha_bootci <- (1 + c(-conf_bootci, conf_bootci))/2
power_val    <- 0.8

# read precomputed raw results data 
fnames_all   <- list.files(res_fdir_raw, full.names = TRUE)
dat <- do.call("rbind", lapply(fnames_all, readRDS))
# fix for a name 
dat <- mutate(dat, name = replace(name, name == "run_result", "run_result_GEE"))

dim(dat)
head(dat)
str(dat)

range(dat$N1)
# [1]  10 250

table(dat$name); Sys.time()
# bootstrap_power_GEE      run_result_GEE   upstrap_power_GEE 
# 31089               31089               31089 
# [1] "2021-06-10 21:00:11 EDT"

# how many out of R=1000 has processed
paste0("arrayjob_idx count = ", length(unique(dat$arrayjob_idx))); Sys.time()
# [1] "arrayjob_idx count = 130"
# [1] "2021-06-10 21:02:15 EDT"


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PART 1: bootstrap to get power estimates (mean, median) and corresponding CI

message("Starting PART 1...")

out_df_all <- data.frame()

# ------------------------------------------------------------------------------
# PART 1A: upstrap_bootstrap_power_diff

name_tmp <- "upstrap_bootstrap_power_diff"
dat_ups   <- dat %>% filter(name == "upstrap_power_GEE")
dat_bots  <- dat %>% filter(name == "bootstrap_power_GEE")

dat_value <- 
  dat_ups %>% 
  left_join(dat_bots, by = c("N1", "arrayjob_idx")) %>%
  mutate(value = value.x - value.y,
         name = name_tmp) %>%
  select(N0 = N0.x, N1 = N1, arrayjob_idx, name, value)

# no need to iterate over N0 
N0_tmp <- 41
dat_i <- 
  dat_value %>% 
  select(arrayjob_idx, N1, value) %>% 
  pivot_wider(names_from = N1, values_from = value) %>%
  select(-arrayjob_idx) %>% 
  as.matrix() 
dat_i_ncol <- ncol(dat_i)
dat_i_nrow <- nrow(dat_i)

# bootstrap the statistics: median, mean 
dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
set.seed(123)
for (b in 1 : B_bootci){
  if (b %% 100 == 0) message(b)
  # resample the dat_value rows
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
  est_median_bootci_lwr = dat_i_bootci_median_lwr,
  est_median_bootci_upr = dat_i_bootci_median_upr,
  est_mean = dat_i_mean,
  est_mean_bootci_lwr = dat_i_bootci_mean_lwr,
  est_mean_bootci_upr = dat_i_bootci_mean_upr
) %>% mutate(
  name = rep(name_tmp, dat_i_ncol),
  N0 = rep(N0_tmp, dat_i_ncol),
  N1 = as.numeric(colnames(dat_i)), 
  .before = everything(),
)
out_df_all <- rbind(out_df_all, out_df_i)
rm(out_df_i, dat_i_bootobj_median, dat_i_bootobj_mean, dat_i)


# ------------------------------------------------------------------------------
# PART 1B: upstrap power 

name_tmp <- "upstrap_power_GEE"

dat_value <- 
  dat %>% filter(name == name_tmp) %>%
  select(N0, N1, arrayjob_idx, name, value)

# no need to iterate over N0 
N0_tmp <- unique(dat_value$N0)[1]

dat_i <- 
  dat_value %>% 
  select(arrayjob_idx, N1, value) %>% 
  pivot_wider(names_from = N1, values_from = value) %>%
  select(-arrayjob_idx) %>% 
  as.matrix() 
dat_i_ncol <- ncol(dat_i)
dat_i_nrow <- nrow(dat_i)

# bootstrap the statistics: median, mean 
dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
set.seed(123)
for (b in 1 : B_bootci){
  if (b %% 100 == 0) message(b)
  # resample the dat_value rows
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
  est_median_bootci_lwr = dat_i_bootci_median_lwr,
  est_median_bootci_upr = dat_i_bootci_median_upr,
  est_mean = dat_i_mean,
  est_mean_bootci_lwr = dat_i_bootci_mean_lwr,
  est_mean_bootci_upr = dat_i_bootci_mean_upr
) %>% mutate(
  name = rep(name_tmp, dat_i_ncol),
  N0 = rep(N0_tmp, dat_i_ncol),
  N1 = as.numeric(colnames(dat_i)), 
  .before = everything(),
)

out_df_all <- rbind(out_df_all, out_df_i)
rm(out_df_i, dat_i_bootobj_median, dat_i_bootobj_mean, dat_i)


# ------------------------------------------------------------------------------
# PART 1C: bootstrap power , run result

for (name_tmp in c("bootstrap_power_GEE", "run_result_GEE")){
  # name_tmp <- "bootstrap_power_GEE"
  dat_value <- 
    dat %>% filter(name == name_tmp) %>%
    select(N0, N1, arrayjob_idx, name, value)
  
  # no need to iterate over N0 
  dat_i <- 
    dat_value %>% 
    select(arrayjob_idx, N1, value) %>% 
    pivot_wider(names_from = N1, values_from = value) %>%
    select(-arrayjob_idx) %>% 
    as.matrix() 
  dat_i_ncol <- ncol(dat_i)
  dat_i_nrow <- nrow(dat_i)
  
  # bootstrap the statistics: median, mean 
  dat_i_bootobj_median <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
  dat_i_bootobj_mean   <- matrix(NA, nrow = B_bootci, ncol = ncol(dat_i))
  set.seed(123)
  for (b in 1 : B_bootci){
    if (b %% 100 == 0) message(b)
    # resample the dat_value rows
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
    est_median_bootci_lwr = dat_i_bootci_median_lwr,
    est_median_bootci_upr = dat_i_bootci_median_upr,
    est_mean = dat_i_mean,
    est_mean_bootci_lwr = dat_i_bootci_mean_lwr,
    est_mean_bootci_upr = dat_i_bootci_mean_upr
  ) %>% mutate(
    name = rep(name_tmp, dat_i_ncol),
    N0 = as.numeric(colnames(dat_i)),
    N1 = as.numeric(colnames(dat_i)), 
    .before = everything(),
  )

  out_df_all <- rbind(out_df_all, out_df_i)
  rm(out_df_i, dat_i_bootobj_median, dat_i_bootobj_mean, dat_i)
}


# ------------------------------------------------------------------------------
# save the aggregated data 

saveRDS(out_df_all, paste0(res_fdir_agg, "/glmm_trt_power_est_bootCI.rds"))
message("Saved the results.")

