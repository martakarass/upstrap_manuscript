#!/usr/bin/env Rscript

#' Notes: 
#' cd $ups 
#' git pull
#' cd $ups/other_experiments/R
#' 
#' ls -l -d *JOB_other2_est*
#' rm JOB_other2_est*
#' 
#' Rnosave 2021-08-06-estimate_power.R -t 1-397 -tc 110 -N JOB_other2_est
#' 
#' rm $ups/other_experiments/results_CL_shared/2021-08-06-estimate_power/*

arg_str <- as.character(Sys.getenv("SGE_TASK_ID"))
arrayjob_idx <- as.numeric(arg_str)

# rm(list = ls()); arrayjob_idx <- 1
library(tidyverse)
message(paste0("arrayjob_idx: ", arrayjob_idx))

out_dir <- paste0(here::here(), "/other_experiments/results_CL_shared/2021-08-06-estimate_power")
dir.create(out_dir)

# parameters
B_boot <- 1000
R_rep  <- 1000 

# grid of parameters
err_sd_grid      <- c(0.01, 0.1, 1, 10)
effsize_tru_grid <- c(0.01)
effsize_tar_grid <- c(0.01, NA) # NA stands for observed
N_tar_grid       <- seq(2000, 16000, by = 2000)
N_obs_grid       <- c(1000, 2000, 4000)

# upstrap -- data frame of grid of parameters
params_df_ups <- list(
  err_sd = err_sd_grid,
  effsize_tru = effsize_tru_grid,
  effsize_tar = effsize_tar_grid,
  N_tar = N_tar_grid,
  N_obs = N_obs_grid
  ) %>% expand.grid() %>%
  mutate(name = "upstrap")

# power.t.test() -- data frame of grid of parameters
params_df_ptt <- list(
  err_sd = err_sd_grid,
  effsize_tru = effsize_tru_grid,
  effsize_tar = effsize_tar_grid,
  N_tar = N_tar_grid,
  N_obs = N_obs_grid
) %>% expand.grid() %>%
  mutate(name = "powerttest")

# "true power" via proportion of rejected null -- data frame of grid of parameters
params_df_tru <- list(
  err_sd = err_sd_grid,
  effsize_tru = effsize_tru_grid,
  effsize_tar = NA,
  N_tar = NA,
  N_obs = N_obs_grid
) %>% expand.grid() %>%
  mutate(name = "true")

# all -- data frame of grid of parameters
params_df <- 
  params_df_ups %>% 
  rbind(params_df_ptt) %>%
  rbind(params_df_tru)
# dim(params_df)
# 396   6


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# util functions

t1 <- Sys.time()

# function to estimate power with upstrap
est_power_upstrap <- function(err_sd, effsize_tru, effsize_tar, N_tar, N_obs, seed_val = 1){
  set.seed(seed_val)
  out_value <- numeric(R_rep)
  for (rr in 1 : R_rep){
    if (rr %% 100 == 0) message(paste0("rr: ", rr, " [", round(rr / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
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

# function to estimate "true power" via proportion of rejections 
est_power_true <- function(err_sd, effsize_tru, effsize_tar, N_tar, N_obs, seed_val = 1){
  set.seed(seed_val)
  out_value <- numeric(R_rep)
  for (rr in 1 : R_rep){
    # print(rr)
    x_rr <- rnorm(n = N_obs, mean = effsize_tru, sd = err_sd)
    out_value[rr] <- (t.test(x_rr)$p.value < 0.05)
  }
  # data frame with the results
  res_df <- data.frame(value = out_value)
  res_df$err_sd = err_sd
  res_df$name = "true"
  res_df$effsize_tru = effsize_tru
  res_df$effsize_tar = NA
  res_df$N_tar = N_tar
  res_df$N_obs = NA
  return(res_df)
}

# function to estimate "true power" via proportion of rejections 
est_power_powerttest <- function(err_sd, effsize_tru, effsize_tar, N_tar, N_obs, seed_val = 1){
  set.seed(seed_val)
  out_value <- numeric(R_rep)
  for (rr in 1 : R_rep){
    if (rr %% 100 == 0) message(paste0("rr: ", rr, " [", round(rr / B_boot * 100, 2), "%], ",  round(as.numeric(Sys.time() - t1, unit = "mins")), " mins"))
    x_rr <- rnorm(n = N_obs, mean = effsize_tru, sd = err_sd)
    if (!is.na(effsize_tar)){ 
      out <- power.t.test(delta = effsize_tar, sd=sd(x_rr), n=N_tar, type="one.sample")
    } else { # observed power
      out <- power.t.test(delta = mean(x_rr), sd=sd(x_rr), n=N_tar, type="one.sample")
    }
    out_value[rr] <- out$power
  }
  # data frame with the results
  res_df <- data.frame(value = out_value)
  res_df$err_sd = err_sd
  res_df$name = "powerttest"
  res_df$effsize_tru = effsize_tru
  res_df$effsize_tar = effsize_tar
  res_df$N_tar = N_tar
  res_df$N_obs = N_obs
  return(res_df)
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# run estimation for current cluster job index 

# pull parameters specific to this job iteration
err_sd      <- params_df[arrayjob_idx, "err_sd"]
effsize_tru <- params_df[arrayjob_idx, "effsize_tru"]
effsize_tar <- params_df[arrayjob_idx, "effsize_tar"]
N_tar       <- params_df[arrayjob_idx, "N_tar"]
N_obs       <- params_df[arrayjob_idx, "N_obs"]
name_tmp    <- params_df[arrayjob_idx, "name"]

# run power estimation 
if (name_tmp == "upstrap"){
  message("Running est_power_upstrap()")
  out <- est_power_upstrap(err_sd, effsize_tru, effsize_tar, N_tar, N_obs)
} else if (name_tmp == "true") {
  message("Running est_power_true()")
  out <- est_power_true(err_sd, effsize_tru, effsize_tar, N_tar, N_obs)
} else  {
  message("Running est_power_powerttest()")
  out <- est_power_powerttest(err_sd, effsize_tru, effsize_tar, N_tar, N_obs)
}

message("EXPERIMENT COMPLETED")

# save results to file 
saveRDS(object = out, file =  paste0(out_dir, "/arrayjob_", arrayjob_idx, ".rds"))
message("DATA SAVED")





