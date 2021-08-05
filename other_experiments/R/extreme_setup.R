#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' cd $ups/other_experiments/R
#' 
#' ls -l -d *adhoc*
#' rm JOB_adhoc_glmm*
#' 
#' Rnosave extreme_setup.R -t 1-50 -tc 50 -N JOB_adhoc_glmm

rm(list = ls())

B_boot <- 1000
R_rep <- 100

err_sd_grid      <- c(0.01, 0.1, 1, 10)
effsize_tru_grid <- c(0.01)
effsize_tar_grid <- c(0.01, NA)
N_tar_grid       <- seq(2000, 16000, by = 2000)
N_obs_grid       <- c(1000, 2000, 4000)
name_grid        <- c("upstrap", "true")

# function to estimate power with upstrap
est_power_upstrap <- function(err_sd, effsize_tru, effsize_tar, N_tar, N_obs, seed_val = 1){
  set.seed(seed_val)
  out_value <- numeric(R_rep)
  for (rr in 1 : R_rep){
    print(rr)
    # simulate sample
    x_rr <- rnorm(n = N_obs, mean = effsize_tru, sd = err_sd)
    # update sample for the target effect size; else, "observed power" is used
    if (!is.na(effsize_tar)){
      x_rr <- x_rr + (effsize_tar - mean(x_rr)) * 1
    }
    out <- rep(NA, B_boot)
    for (bb in 1 : B_boot){
      x_rr_bb <- sample(x_rr, size = N_tar, replace = TRUE)
      out[bb] <- (t.test(x_rr_bb)$p.value < 0.05)
    }
    out_value[rr] <- mean(out)
  }
  # data frame with the results
  res_df <- data.frame(value = out_value)
  res_df$err_sd = err_sd
  res_df$name = "upstrap"
  res_df$effsize_tru = effsize_tru
  res_df$effsize_tar = effsize_tar
  res_df$N_tar = N_tar
  res_df$N_obs = N_obs
  return(res_df)
}

err_sd      <- 10
effsize_tru <- 0.01
effsize_tar <- NA
N_tar       <- 2000
N_obs       <- 1000
out <- est_power_upstrap(err_sd, effsize_tru, effsize_tar, N_tar, N_obs)

# 1 -> NA 
summary(out$value)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0460  0.0835  0.1930  0.2727  0.3595  0.9910 

# 1 -> 0.01 
summary(out$value)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0480  0.0670  0.0730  0.0727  0.0790  0.0920 


# 10 -> NA 
summary(out$value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04000 0.08775 0.13250 0.24172 0.34050 0.94700 

# 10 -> 0.01 
summary(out$value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03600 0.04575 0.04900 0.04980 0.05400 0.06500 













