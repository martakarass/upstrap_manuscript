
#' This script asks to answer a few questions on how standard deviation estimator,
#' s, is estimating the true sigma comparing these methods: 
#' 
#' (1) basically using sd() on observed sample of size N1 
#' (2) bootstrapping observed sample of size N1, using mean(sd(boot_resample_i)) 
#' (3) upstrapping observed sample of size N1 to size N0 (take head) to size N1,
#'     then using mean(sd(upboot_resample_i)) 
 
rm(list = ls())
library(here)
library(tidyverse)
options(dplyr.summarise.inform = FALSE) 
library(matrixStats)
library(data.table)

# dir to save results 
res_fdir <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-01-05-sd_estimation_test")

# params
sigma_grid <- c(0.5, 1, 2, 5, 10)
N0_grid    <- c(30, 50, 100, 200)
N1_grid    <- c(30, 50, 100, 200)
N1_max     <- max(N1_grid)
rep_n  <- 1000 * 100
B_boot <- 1000

sigma_grid_l <- length(sigma_grid)
N1_grid_l    <- length(N1_grid)
N0_grid_l    <- length(N0_grid)


# ------------------------------------------------------------------------------
# function to compute cumulative var
cumvar <- function (x) {
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  return(v)
}

cumsd <- function (x) {
  v <- cumvar(x)
  v <- sqrt(v)
  return(v)
}


# ------------------------------------------------------------------------------
# (1)

set.seed(123)

out_df_agg_ALL <- data.frame()

t1 <- Sys.time()

for (sigma_idx in 1 : sigma_grid_l){
  for (N0_grid_idx in 1 : N0_grid_l){ # sigma_idx <- 1; N0_grid_idx <- 1

    sigma_TMP <- sigma_grid[sigma_idx]
    N0_TMP    <- N0_grid[N0_grid_idx]
    
    message(paste0("STARTING:  sigma_TMP: ", sigma_TMP, ", N0_TMP: ", N0_TMP))

    # iterate over number of repetitions
    out_dt <- data.table()
    for (rep_i in 1 : rep_n){
      # sample current "observed sample" 
      sample_i <- rnorm(n = N0_TMP, mean = 0, sd = sigma_TMP)
      # (1) basically using sd() on observed sample of size N0 
      est_sd_sample_i <- sd(sample_i)
      # (2) bootstrapping observed sample of size N0, using mean(sd(boot_resample_i)) 
      boot_resamples_i <- matrix(sample(x = sample_i, size = (B_boot * N1_max), replace = TRUE), nrow = B_boot)
      boot_cumsd <- apply(boot_resamples_i, 1, cumsd) # apply's output is [N1_max x B_boot]
      est_sd_sample_i_N1 <- matrixStats::rowMeans2(boot_cumsd, na.rm = TRUE)[N1_grid]
      # append current repetition results
      out_dt_i <- data.table(
        N0 = rep(N0_TMP, 1 + N1_grid_l),
        N1 = c(N0_TMP, N1_grid),
        sigma_true = rep(sigma_TMP, 1 + N1_grid_l),
        rep_i = rep(rep_i, 1 + N1_grid_l),
        group = c("sample_sd", rep("boot_samples_sd", N1_grid_l)),
        sd_est = c(est_sd_sample_i, est_sd_sample_i_N1)
      )
      out_dt = rbindlist(list(out_dt, out_dt_i))
    }
    
    # aggregate and append (sigma, N0)-specific results 
    out_df_agg_TMP <- 
      out_dt %>% as.data.frame() %>%
      group_by(N0, N1, sigma_true, group) %>%
      summarise(sd_est_sd = sd(sd_est),
                sd_est_mean = mean(sd_est),
                sd_est_median = median(sd_est),
                rep_n = n()) %>%
      as.data.frame()
    
    out_df_agg_ALL <- base::rbind(out_df_agg_ALL, out_df_agg_TMP)
    
    t_diff <- round(as.numeric(Sys.time() - t1, unit = "secs") / (60 * 60), 5)
    message(paste0("COMPLETED: sigma_TMP: ", sigma_TMP, ", N0_TMP: ", N0_TMP, " | hours elapsed: ", t_diff))
  }
}

out_df_agg_ALL <- as.data.frame(out_df_agg_ALL)

t2 <- Sys.time()
t2-t1

# Time difference of 47.70957 secs

out_df_fpath <- paste0(res_fdir, "/res_repn_", rep_n, "_boot_", B_boot, "_sd_estimation_test")
saveRDS(out_df_agg_ALL, out_df_fpath)

# plt_df <- out_df_agg_ALL %>% filter(N0 == N1) %>%
#   pivot_longer(cols = c(sd_est_median, sd_est_mean)) %>%
#   mutate(group2 = paste0(group, "_", name))
# plt_df_hline <- plt_df %>% select(sigma_true) %>% distinct() %>% mutate(yintercept = sigma_true)
# 
# ggplot(plt_df, aes(x = N0, y = value / sigma_true, color = group2)) +
#   geom_hline(yintercept = 1) +
#   facet_wrap(~ sigma_true, ncol = 1) +
#   geom_point()




