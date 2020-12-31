
#' This script asks to answer a few questions on how power estimation behaves 
#' when: 
#' 
#' (1) we know theoretical result
#' (2) we simulate independent samples (nested independent, or independent-independent)
#' (3) we use sample estimate   

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)

# dir to save results 
res_fdir <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2020-12-31-twosample_ttest")

# params
gg_base_size <- 12
N0 <- 50
N1_min <- 50
N1_max <- 250 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)
N0_grid <- c(50, 100, 150, 200)

# data generating model
mu <- 0.3
simga2 <- 1

# number of repetitions
rep_n <- 100
B_boot <- 1000


# ------------------------------------------------------------------------------
# (1)

out <- sapply(N1_grid, function(n_tmp) power.t.test(n = n_tmp, delta = mu, sd = 1, type = "two.sample")$power)
out_df <- data.frame(N1 = N1_grid, power_est = out)

# ggplot(out_df, aes(x = N1, y = power_est)) + 
#   geom_line() + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
#   scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
#   geom_hline(yintercept = 0.8, color = "blue") + 
#   theme_grey(base_size = gg_base_size) + 
#   labs(title = "two-sample t-test: theoretical power") 

out_df_fpath <- paste0(res_fdir, "/res_theoret")
saveRDS(out_df, out_df_fpath)


# ------------------------------------------------------------------------------
# (2)

t1 <- Sys.time()

# (a) we simulate independent samples (nested independent)
set.seed(123)
mat_out_2a <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
x_mat_gr1  <- matrix(rnorm(n = rep_n * N1_max, mean = 0), ncol = N1_max)
x_mat_gr2  <- matrix(rnorm(n = rep_n * N1_max, mean = mu), ncol = N1_max)
for (i in 1 : rep_n){ # i <- 1; j <- 10
  for (j in 1 : N1_grid_l){
    N1_tmp <- N1_grid[j]
    x_TMP_gr1 <- x_mat_gr1[i, 1 : N1_tmp]
    x_TMP_gr2 <- x_mat_gr2[i, 1 : N1_tmp]
    mat_out_2a[i, j] <- t.test(x_TMP_gr1, x_TMP_gr2, var.equal = TRUE)$p.value < 0.05
  }
}
out_vec_2a <- colMeans(mat_out_2a)
message("Completed out_vec_2a")

# (b) we simulate independent samples (independent-independent)
set.seed(123)
mat_out_2b <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
for (i in 1 : rep_n){ # i <- 1; j <- 10
  for (j in 1 : N1_grid_l){
    N1_tmp <- N1_grid[j]
    x_TMP_gr1 <- rnorm(n = N1_tmp, mean = 0)
    x_TMP_gr2 <- rnorm(n = N1_tmp, mean = mu)
    mat_out_2b[i, j] <- t.test(x_TMP_gr1, x_TMP_gr2, var.equal = TRUE)$p.value < 0.05
  }
}
out_vec_2b <- colMeans(mat_out_2b)
message("Completed out_vec_2bs")

out_df_2 <- 
  data.frame(N1 = N1_grid, rv_nested = out_vec_2a, rv_unique = out_vec_2b) %>%
  pivot_longer(cols = c(rv_nested, rv_unique), 
               names_to = "group", 
               values_to = "power_est")

# ggplot(out_df_2, aes(x = N1, y = power_est, color = group)) + 
#   geom_line() + 
#   geom_line(data = out_df,  aes(x = N1, y = power_est), color = "black") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
#   scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
#   geom_hline(yintercept = 0.8, color = "blue") + 
#   theme_grey(base_size = gg_base_size) + 
#   labs(title = "two-sample t-test: empirical power -- gold standard") 

out_df_fpath <- paste0(res_fdir, "/res_repn_", rep_n, "_boot_", B_boot, "_indepsampling")
saveRDS(out_df_2, out_df_fpath)



# ------------------------------------------------------------------------------
# (3)

# compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
cum_rejectH0_twosample_ttest <- function(vals1, vals2){
  vals_cumnx      <- seq_along(vals1)
  vals_cumdf      <- 2 * vals_cumnx - 2
  vals1_cummean    <- cumsum(vals1) / vals_cumnx
  vals2_cummean    <- cumsum(vals2) / vals_cumnx
  vals_cummeandiff <- vals1_cummean - vals2_cummean
  vals1_sumsquares <- (cumsum(vals1 ^ 2) - cumsum(vals1) ^ 2 / seq_along(vals1))
  vals2_sumsquares <- (cumsum(vals2 ^ 2) - cumsum(vals2) ^ 2 / seq_along(vals2))
  vals_cumvarpool  <- (vals1_sumsquares + vals2_sumsquares) / vals_cumdf
  vals_cumsizeinvsum  <- (1 / vals_cumnx) + (1 / vals_cumnx)
  vals_cumstderr   <- sqrt(vals_cumvarpool * vals_cumsizeinvsum)
  vals_cumststat   <- vals_cummeandiff / vals_cumstderr
  vals_cumpval     <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
  vals_reject_H0   <- vals_cumpval < 0.05
  return(vals_reject_H0)
}

N0_vec        <- numeric()
N1_vec        <- numeric()
group_vec     <- numeric()
power_est_vec <- numeric()

set.seed(123)
for (N0 in N0_grid){ # N0 <- 50; i <- 1
  message(N0)
  N1_grid <- N0 : N1_max
  N1_grid_l <- length(N1_grid)
  # set of samples for current N0
  x_mat_gr1  <- matrix(rnorm(n = rep_n * N0, mean = 0), nrow = rep_n)
  x_mat_gr2  <- matrix(rnorm(n = rep_n * N0, mean = mu), nrow = rep_n)
  mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  for (i in 1:rep_n){
    sample_i_gr1 <- x_mat_gr1[i, ]
    sample_i_gr2 <- x_mat_gr2[i, ]
    # sample_i: generate power estimate with power.t.test()
    sample_i_meandiff   <-  mean(sample_i_gr1) - mean(sample_i_gr2)
    sample_i_var_pooled <- (var(sample_i_gr1) * (N0 - 1) + var(sample_i_gr2) * (N0 - 1)) / (N0 + N0 - 2)
    sample_i_sd_pooled  <- sqrt(sample_i_var_pooled)
    mat_out_powerttest[i, ] <- sapply(N1_grid, function(N1_tmp){
      (power.t.test(n = N1_tmp, 
                    delta = sample_i_meandiff, 
                    sd = sample_i_sd_pooled, 
                    sig.level = 0.05, type = "two.sample", alternative = "two.sided"))$power
    })
    # sample_i: generate power estimate with upstrap()
    boot_resamples_i_gr1 <- matrix(sample(x = sample_i_gr1, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_gr2 <- matrix(sample(x = sample_i_gr2, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_rejectH0 <- lapply(1:B_boot, function(k) cum_rejectH0_twosample_ttest(
      boot_resamples_i_gr1[k, ],
      boot_resamples_i_gr2[k, ]
    ))
    boot_resamples_i_rejectH0 <- do.call(rbind, boot_resamples_i_rejectH0)
    mat_out_upstrap[i, ] <- apply(boot_resamples_i_rejectH0[, N1_grid], 2, mean, na.rm = TRUE)
  }
  out_powerttest_mean   <- apply(mat_out_powerttest, 2, mean)
  out_powerttest_median <- apply(mat_out_powerttest, 2, median)
  out_upstrap_mean      <- apply(mat_out_upstrap, 2, mean)
  out_upstrap_median    <- apply(mat_out_upstrap, 2, median)
  
  N0_vec        <- c(N0_vec, rep(N0, 4 * N1_grid_l))
  N1_vec        <- c(N1_vec, rep(N1_grid, 4))
  group_vec     <- c(group_vec, 
                     rep("powerttest_mean", N1_grid_l), 
                     rep("powerttest_median", N1_grid_l),
                     rep("upstrap_mean", N1_grid_l), 
                     rep("upstrap_median", N1_grid_l))
  power_est_vec <- c(power_est_vec, 
                     out_powerttest_mean,
                     out_powerttest_median,
                     out_upstrap_mean,
                     out_upstrap_median)
}


out_df_3 <- 
  data.frame(N0 = N0_vec, N1 = N1_vec, group = group_vec, power_est = power_est_vec) 

out_df_fpath <- paste0(res_fdir, "/res_repn_", rep_n, "_boot_", B_boot, "_powerttest_upstrap")
saveRDS(out_df_3, out_df_fpath)
message("Completed out_df_3")

# ggplot(out_df_3, aes(x = N1, y = power_est, color = group)) + 
#   geom_line() + 
#   geom_line(data = out_df,  aes(x = N1, y = power_est), color = "black") + 
#   facet_grid(N0 ~ .) + 
#   # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
#   scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
#   geom_hline(yintercept = 0.8, color = "blue") + 
#   theme_grey(base_size = gg_base_size) + 
#   labs(title = "two-sample t-test: powerttest/upstrap-estimated power") 

t2 <- Sys.time()
message(t2 - t1)

# Time difference of 7.33839 mins



