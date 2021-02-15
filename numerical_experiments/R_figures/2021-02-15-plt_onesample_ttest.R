
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
res_fdir <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2020-12-30-onesample_ttest")

# experiment parameteres
gg_base_size <- 12
N0_grid <- c(30, 50, 100, 150)
N1_min <- 30
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu     <- 0.3
simga2 <- 1

# number of repetitions pf experiment 
rep_n   <- 10000
# number of boot repetitions within one experiment, one setup
B_boot  <- 1000


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# (1) get the estimates using theoretical (gold standard) results

out <- sapply(N1_grid, function(n_tmp) power.t.test(n = n_tmp, delta = mu, sd = 1, type = "one.sample")$power)
out_df <- data.frame(N1 = N1_grid, power_est = out)
ggplot(out_df, aes(x = N1, y = power_est)) + 
  geom_line() + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
  geom_hline(yintercept = 0.8, color = "blue") + 
  theme_grey(base_size = gg_base_size) + 
  labs(title = "one-sample t-test: theoretical power") 

out_df_fpath <- paste0(res_fdir, "/res_theoret")
saveRDS(out_df, out_df_fpath)


# ------------------------------------------------------------------------------
# (2) get the estimates using sampling of independent samples

t1 <- Sys.time()

# (a) we simulate independent samples (nested independent)
set.seed(123)
x_mat <- matrix(rnorm(n = rep_n * N1_max, mean = mu), ncol = N1_max)
out_vec_2a <- numeric()
for (i in 1 : N1_grid_l){
  n_tmp <- N1_grid[i]
  out_vec_2a[i] <- mean(apply(x_mat[, 1 : n_tmp], 1, function(row_i) t.test(row_i)$p.value < 0.05))
}

# (b) we simulate independent samples (independent-independent)
set.seed(123)
out_vec_2b <- numeric()
for (i in 1 : N1_grid_l){
  n_tmp <- N1_grid[i]
  x_mat <- matrix(rnorm(n = rep_n * n_tmp, mean = mu), ncol = n_tmp)
  out_vec_2b[i] <- mean(apply(x_mat, 1, function(row_i) t.test(row_i)$p.value < 0.05))
}

out_df_2 <- 
  data.frame(N1 = N1_grid, rv_nested = out_vec_2a, rv_unique = out_vec_2b) %>%
  pivot_longer(cols = c(rv_nested, rv_unique), 
               names_to = "group", 
               values_to = "power_est")
ggplot(out_df_2, aes(x = N1, y = power_est, color = group)) + 
  geom_line() + 
  geom_line(data = out_df,  aes(x = N1, y = power_est), color = "black") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
  geom_hline(yintercept = 0.8, color = "blue") + 
  theme_grey(base_size = gg_base_size) + 
  labs(title = "one-sample t-test: theoretical power") 

out_df_fpath <- paste0(res_fdir, "/res_repn_", rep_n, "_boot_", B_boot, "_indepsampling")
out_df_fpath
saveRDS(out_df_2, out_df_fpath)


# ------------------------------------------------------------------------------
# (2) get the estimates using (a) power.t.test(), (b) usptrap

# function to compute cumulative var
cumvar <- function (x, sd = FALSE) {
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  # if (sd) v <- sqrt(v)
  v
}

# compute "cumulative" t.test results (1 - reject H0, 0 -- do not reject H0)  
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

N0_vec        <- numeric()
N1_vec        <- numeric()
group_vec     <- numeric()
power_est_vec <- numeric()

set.seed(123)
for (N0 in N0_grid){ # N0 <- 50; i <- 1
  print(N0)
  N1_grid <- N0 : N1_max
  N1_grid_l <- length(N1_grid)
  # set of samples for current N0
  x_mat <- matrix(rnorm(n = rep_n * N0, mean = mu), nrow = rep_n)
  mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = N1_grid_l)
  for (i in 1:rep_n){
    sample_i <- x_mat[i, ]
    # sample_i: generate power estimate with power.t.test()
    sample_i_mean <- mean(sample_i)
    sample_i_sd   <- sd(sample_i)
    mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
      (power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
    })
    # sample_i: generate power estimate with upstrap()
    boot_resamples_i <- matrix(sample(x = sample_i, size = (B_boot * N1_max), replace = TRUE), 
                               nrow = B_boot, ncol = N1_max)
    boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
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

# save the results 
out_df_3 <- 
  data.frame(N0 = N0_vec, N1 = N1_vec, group = group_vec, power_est = power_est_vec) 

out_df_fpath <- paste0(res_fdir, "/res_repn_", rep_n, "_boot_", B_boot, "_powerttest_upstrap")
saveRDS(out_df_3, out_df_fpath)

# ggplot(out_df_3, aes(x = N1, y = power_est, color = group)) + 
#   geom_line() + 
#   geom_line(data = out_df,  aes(x = N1, y = power_est), color = "black") + 
#   facet_grid(N0 ~ .) + 
#   # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
#   scale_x_continuous(limits = c(N1_min, N1_max), breaks = seq(N1_min, N1_max, by = 10)) + 
#   geom_hline(yintercept = 0.8, color = "blue") + 
#   theme_grey(base_size = gg_base_size) + 
#   labs(title = "one-sample t-test: theoretical power") 

t2 <- Sys.time()
t2-t1




