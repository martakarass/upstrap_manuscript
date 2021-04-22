
# ------------------------------------------------------------------------------
# PLOT 1 

rm(list = ls())

# simulation parameters
N      <- 30
M_min  <- 30
M_max  <- 200
M_grid <- M_min : M_max
M_grid_l <- length(M_grid)
R <- 1000
B <- 1000
mu <- 0.3
sigma2 <- 1

# objects to store simulation results
rep_idx_vec       <- numeric()
M_sample_size_vec <- numeric()
powerttest_vec    <- numeric()
upstrap_vec       <- numeric()
ttest_vec         <- numeric()

# function to compute "cumulative" t.test result (1 - reject H0, 0 -- do not reject H0)  
vals_cum_reject_H0 <- function(vals){
  vals_cumnx      <- seq_along(vals)
  vals_cumdf      <- vals_cumnx - 1
  vals_cummean    <- cumsum(vals) / vals_cumnx
  vals_cumvar     <- (cumsum(vals ^ 2) - cumsum(vals) ^ 2 / vals_cumnx) / (vals_cumnx - 1)
  vals_cumstderr  <- sqrt(vals_cumvar / vals_cumnx)
  vals_cumststat  <- vals_cummean / vals_cumstderr
  vals_cumpval    <- 2 * pt(-abs(vals_cumststat), vals_cumdf)
  vals_reject_H0  <- vals_cumpval < 0.05
  return(vals_reject_H0)
}

for (rep_idx in 1:R){
  print(rep_idx)
  set.seed(rep_idx)
  # simulate sample of full size (M max) and its "observed" subset 
  x_M_max <- rnorm(M_max, mean = mu, sd = sigma2)
  x_N <- x_M_max[1:N]
  # test power: power.t.test
  out <- power.t.test(n = M_grid, delta = mean(x_N), sd = sd(x_N), type = "one.sample")$power
  powerttest_vec <- c(powerttest_vec, out); rm(out)
  # test power: upstrap
  x_N_resamples <- matrix(sample(x_N, size = M_max * B, replace = TRUE), ncol = M_max, nrow = B)
  x_N_resamples_rejectH0 <- t(apply(x_N_resamples, 1, vals_cum_reject_H0))
  out <- apply(x_N_resamples_rejectH0[, M_grid], 2, mean, na.rm = TRUE)
  upstrap_vec <- c(upstrap_vec, out); rm(out)
  # test result: t.test 
  out <- as.numeric(vals_cum_reject_H0(x_M_max)[M_grid])
  ttest_vec <- c(ttest_vec, out); rm(out)
  # meta data 
  rep_idx_vec <- c(rep_idx_vec, rep(rep_idx, M_grid_l))
  M_sample_size_vec <- c(M_sample_size_vec, M_grid)
}

# combine results into a data frame
out_df <- data.frame(
  rep_idx = rep_idx_vec,
  M_sample_size = M_sample_size_vec,
  powerttest = powerttest_vec,
  upstrap = upstrap_vec,
  ttest = ttest_vec
)

# dir to save results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
saveRDS(out_df, file = paste0(results_dir, "/out_df.rds"))


# ------------------------------------------------------------------------------
# PLOT 2

rm(list = ls())
source(here::here("numerical_experiments/R/config_utils.R"))

# simulation parameters
N      <- 50
M_min  <- 30
M_max  <- 200
M_grid <- M_min : M_max
M_grid_l <- length(M_grid)
R <- 1000
# B <- 1000
mu <- 0.3
sigma2 <- 1

# simulate draws from sampling distribution of sample mean 
set.seed(1)
xbar_obs_vec   <- rnorm(R, mean = mu, sd = (sqrt(sigma2)/sqrt(N)))
xbar_obs_vec_l <- length(xbar_obs_vec)

# simulate draws from sampling distribution of sample standard deviation (s) 
set.seed(1)
s_obs_vec  <- sqrt(rchisq(R, df = (N - 1), ncp = 0)) * sqrt(1 / (N-1)) * 1
# s_obs_vec   <- rnorm(R, mean = s_mean, sd = s_sd)
s_obs_vec_l <- length(s_obs_vec)

# objects to store simulation results
rep_idx_vec       <- numeric()
M_sample_size_vec <- numeric()
powerttest_vec    <- numeric()
example_vec       <- numeric()

for (rep_idx in 1:R){ # rep_idx <- 1 
  print(rep_idx)
  xbar_obs_i <- xbar_obs_vec[rep_idx]
  s_obs_i    <- s_obs_vec[rep_idx]
  # get power vals
  out_a <- power.t.test(n = M_grid, delta = xbar_obs_i, sd = sqrt(sigma2), type = "one.sample")$power
  out_b <- power.t.test(n = M_grid, delta = mu,         sd = s_obs_i, type = "one.sample")$power
  out_c <- power.t.test(n = M_grid, delta = xbar_obs_i, sd = s_obs_i, type = "one.sample")$power
  # append
  rep_idx_vec       <- c(rep_idx_vec, rep(rep_idx, M_grid_l * 3))
  M_sample_size_vec <- c(M_sample_size_vec, M_grid, M_grid, M_grid)
  powerttest_vec    <- c(powerttest_vec, out_a, out_b, out_c)
  example_vec       <- c(example_vec, rep("ex_a", M_grid_l), rep("ex_b", M_grid_l), rep("ex_c", M_grid_l))
} 

# combine results into a data frame
out_df <- data.frame(
  rep_idx = rep_idx_vec,
  M_sample_size = M_sample_size_vec,
  powerttest = powerttest_vec,
  example = example_vec
)

# dir to save results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
saveRDS(out_df, file = paste0(results_dir, "/out_df_PLOT2.rds"))



# ------------------------------------------------------------------------------
# PLOT 3

rm(list = ls())
source(here::here("numerical_experiments/R/config_utils.R"))
library(chi)

# simulation parameters
N   <- 50
R   <- 100000
mu  <- 0.3
sigma2 <- 1
M_k <- 100

xbar_obs_vec <- seq(-0.1, 0.7, length.out = 1000)
s_obs_vec <- seq(0.75, 1.25, length.out = 1000)

# # simulate draws from sampling distribution of sample mean 
# set.seed(1)
# xbar_obs_vec   <- rnorm(R, mean = mu, sd = (sqrt(sigma2)/sqrt(N)))
# 
# # simulate draws from sampling distribution of sample standard deviation (s) 
# set.seed(1)
# s_obs_vec  <- sqrt(rchisq(R, df = (N - 1), ncp = 0)) * sqrt(1 / (N-1)) * 1

# objects to store simulation results
arg_val  <- c(
  xbar_obs_vec, 
  s_obs_vec)
func_val <- c(
  power.t.test(n = M_k, delta = xbar_obs_vec, sd = sqrt(sigma2), type = "one.sample")$power,
  power.t.test(n = M_k, delta = mu, sd = s_obs_vec, type = "one.sample")$power
)
dens_val <- c(
  power.t.test(n = M_k, delta = xbar_obs_vec, sd = sqrt(sigma2), type = "one.sample")$power,
  power.t.test(n = M_k, delta = mu, sd = s_obs_vec, type = "one.sample")$power
)
arg_name <- c(
  rep("xbar", R),
  rep("s", R)
)
out_df <- data.frame(
  arg_val, 
  func_val,
  dens_val,
  arg_name)
head(out_df)

# dir to save results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
saveRDS(out_df, file = paste0(results_dir, "/out_df_PLOT3.rds"))

