
rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-06-28-onesample_ttest_agg")
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-28-onesample_ttest_raw")
res_fdir_plt  <- paste0(here::here(), "/numerical_experiments/results_figures/2021-06-28-onesample_ttest")
# remove dirs if exist to make a room for new ones
unlink(res_fdir_agg, recursive = TRUE, force = TRUE)
unlink(res_fdir_raw, recursive = TRUE, force = TRUE)
unlink(res_fdir_plt, recursive = TRUE, force = TRUE)
# create dirs if any does not exist
dir.create(path = res_fdir_agg)
dir.create(path = res_fdir_raw)
dir.create(path = res_fdir_plt)
message(paste0("dir.exists(path = res_fdir_agg): ", dir.exists(path = res_fdir_agg)))
message(paste0("dir.exists(path = res_fdir_raw): ", dir.exists(path = res_fdir_raw)))
message(paste0("dir.exists(path = res_fdir_plt): ", dir.exists(path = res_fdir_plt)))

# experiment parameters
N0 <- 45
N1_min <- 5
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu0   <- 0.3
simga2 <- 1

# number of repetitions of experiment 
R_rep   <- 1000 
# number of boot repetitions within one experiment, one setup
B_boot  <- 1000

# effect of interest
mu1_grid <- c(0.2, 0.3, 0.4)


# ------------------------------------------------------------------------------
# HELPER FUNCTIONS

# function to compute cumulative var
cumvar <- function (x, sd = FALSE) {
    n <- seq_along(x)
    v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
    # if (sd) v <- sqrt(v)
    v
}

# function to compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
vals_cum_reject_H0 <- function(vals){
    vals_cumnx      <- seq_along(vals)
    vals_cumdf      <- vals_cumnx - 1
    vals_cummean    <- cumsum(vals) / vals_cumnx
    vals_cumvar     <- cumvar(vals)
    vals_cumstderr  <- sqrt(vals_cumvar / vals_cumnx)
    vals_cumststat  <- vals_cummean / vals_cumstderr
    vals_cumpval    <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
    vals_reject_H0  <- vals_cumpval < 0.05
    return(vals_reject_H0)
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATION

for (arrayjob_idx in 1 : R_rep){
    
    set.seed(arrayjob_idx)
    message(paste0("arrayjob_idx: ", arrayjob_idx))
    
    # simulate sample observed in this experiment repetition
    sample_i <- rnorm(n = N0, mean = mu0)
    sample_i_mean <- mean(sample_i)
    sample_i_sd   <- sd(sample_i)
    
    mat_out_all <- data.frame()
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH power.t.test()
    
    value <- sapply(N1_grid, function(n_tmp){
        out_test <- power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
        out_test$power
    })
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$name          <- "powerttest_power"
    mat_out_tmp$mu0           <- mu0
    mat_out_tmp$mu1           <- NA
    mat_out_tmp$value <- value
    mat_out_all <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
    
    
    # ESTIMATE POWER WITH power.t.test(), set mu1
    for (mu1 in mu1_grid){
        value <- sapply(N1_grid, function(n_tmp){
            out_test <- power.t.test(n = n_tmp, delta = mu1, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
            out_test$power 
        })
        mat_out_tmp               <- data.frame(N1 = N1_grid)
        mat_out_tmp$N0            <- rep(N0, N1_grid_l)
        mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
        mat_out_tmp$name          <- "powerttest_power"
        mat_out_tmp$mu0           <- mu0
        mat_out_tmp$mu1           <- mu1
        mat_out_tmp$value <- value
        mat_out_all <- rbind(mat_out_all, mat_out_tmp)
        rm(mat_out_tmp, value)
    }
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH upstrap
    
    sample_i_updated <- sample_i 
    boot_resamples_i <- matrix(sample(x = sample_i_updated, size = (B_boot * N1_max), replace = TRUE), nrow = B_boot, ncol = N1_max)
    boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
    value <- apply(boot_resamples_i_rejectH0[, N1_grid], 2, mean, na.rm = TRUE)
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$name          <- "upstrap_power"
    mat_out_tmp$mu0           <- mu0
    mat_out_tmp$mu1           <- NA
    mat_out_tmp$value <- value
    mat_out_all <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
    
    
    # ESTIMATE POWER WITH upstrap, set mu1
    for (mu1 in mu1_grid){
        sample_i_updated <- sample_i + (mu1 - sample_i_mean)
        boot_resamples_i <- matrix(sample(x = sample_i_updated, size = (B_boot * N1_max), replace = TRUE), nrow = B_boot, ncol = N1_max)
        boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
        value <- apply(boot_resamples_i_rejectH0[, N1_grid], 2, mean, na.rm = TRUE)
        mat_out_tmp               <- data.frame(N1 = N1_grid)
        mat_out_tmp$N0            <- rep(N0, N1_grid_l)
        mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
        mat_out_tmp$name          <- "upstrap_power"
        mat_out_tmp$mu0           <- mu0
        mat_out_tmp$mu1           <- mu1
        mat_out_tmp$value <- value
        mat_out_all <- rbind(mat_out_all, mat_out_tmp)
        rm(mat_out_tmp, value)
    }
    
    
    # ------------------------------------------------------------------------------
    # SAVE TO FILE 
    out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
    saveRDS(object = mat_out_all, file = out_fpath_raw)
}




