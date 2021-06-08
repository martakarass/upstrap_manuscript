
# ------------------------------------------------------------------------------
# PLOT 1 

rm(list = ls())

# simulation parameters
N      <- 25
M_min  <- 5
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
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-05-29-onesample_ttest_aggegating_comparison")
saveRDS(out_df, file = paste0(results_dir, "/out_df.rds"))

