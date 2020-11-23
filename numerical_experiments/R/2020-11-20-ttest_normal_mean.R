rm(list = ls())

project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
# project_dir <- "/users/mkaras/_PROJECTS/upstrap_manuscript"
library(data.table)

# params 
mu        <- 0.1 
n0        <- 200
rep_n     <- 10000
B_boot    <- 10000
n_grid    <- n0:1200
n_grid_max <- max(n_grid)
print(paste0("n0: ", n0, ", rep_n: ", rep_n, ", B_boot: ", B_boot))

# function to compute cumulative var
cumvar <- function (x, sd = FALSE) {
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  # if (sd) v <- sqrt(v)
  v
}

# compute "cumulative t.test results (1 - reject H0, 0 -- do not reject H0)"  
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

# objects to store simulation resuls
mat_out_powerttest <- matrix(NA, nrow = rep_n, ncol = length(n_grid))
mat_out_upstrap    <- matrix(NA, nrow = rep_n, ncol = length(n_grid))
dim(mat_out_powerttest)
dim(mat_out_upstrap)

mat_out_powerttest_fpath <- paste0(project_dir, "/numerical_experiments/results/2020-11-20-ttest_normal_mean/mat_out_powerttest.csv")
mat_out_upstrap_fpath    <- paste0(project_dir, "/numerical_experiments/results/2020-11-20-ttest_normal_mean/mat_out_upstrap.csv")

t1 <- Sys.time()
set.seed(123)
for (i in 1:rep_n){
  print(i)
  # simulate sample 
  sample_i <- rnorm(n = n0, mean = mu, sd = 1)
  # sample_i: generate power estimate with power.t.test()
  sample_i_mean <- mean(sample_i)
  sample_i_sd   <- sd(sample_i)
  mat_out_powerttest[i, ] <- sapply(n_grid, function(n_tmp){
    (power.t.test(n = n_tmp, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
  })
  mat_out_powerttest[i, ] <- round(mat_out_powerttest[i, ], 5)
  # sample_i: generate power estimate with upstrap()
  # generate the samples
  boot_resamples_i <- matrix(sample(x = sample_i, size = (B_boot * n_grid_max), replace = TRUE), 
                             nrow = B_boot, ncol = n_grid_max)
  boot_resamples_i_rejectH0 <- t(apply(boot_resamples_i, 1, vals_cum_reject_H0))
  mat_out_upstrap[i, ] <- apply(boot_resamples_i_rejectH0[, n_grid], 2, mean)
  mat_out_upstrap[i, ] <- round(mat_out_upstrap[i, ], 5)
  # save to file every 100
  if (i %% 500 == 0){
    fwrite(as.data.table(mat_out_powerttest), mat_out_powerttest_fpath)
    fwrite(as.data.table(mat_out_upstrap),    mat_out_upstrap_fpath)
    secs_passed <- round(as.numeric(Sys.time() - t1, unit = "secs"))
    paste0("i: ", i, " [", round(i / rep_n * 100, 2), "%], ", secs_passed, " secs")
  }
  
}
t2 <- Sys.time()

# final save
fwrite(as.data.table(mat_out_powerttest), mat_out_powerttest_fpath)
fwrite(as.data.table(mat_out_upstrap),    mat_out_upstrap_fpath)
secs_passed <- round(as.numeric(t1 - t2, unit = "secs"))
paste0("i: ", i, " [", round(i / rep_n * 100, 2), "%], ", secs_passed, " secs")





