
#' @description 
#' Script to estimate CI coverage of parameter true value in the experiment
#' of estimating N(0,1) sample mean. 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>

## -----------------------------------------------------------------------------

rm(list = ls())

# params
project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
out_fname <- "2020-09-13-estimate_ci_coverage_N_results.txt"
N_seq <- seq(from = 5, to = 250, by = 5)
R <- 10000


## -----------------------------------------------------------------------------

out_fpath <- paste0(project_dir, "/numerical_experiments/results/", out_fname)

library(data.table)

# testing purpose
get_boot_CI <- function(vals, B_boot = 10000){
  vals_l <- length(vals)
  boot_samples <- matrix(sample(x = vals, size = (B_boot * vals_l), replace = TRUE), nrow = B_boot, ncol = vals_l)
  boot_sample_means <- apply(boot_samples, 1, mean)
  m_boot  <- mean(boot_sample_means) 
  sd_boot <- sd(boot_sample_means) 
  CI_boot_1 <- m_boot + c(-1, 1) * qnorm(1 - 0.025) * sd_boot 
  CI_boot_2 <- m_boot + c(-1, 1) * qt(1 - 0.025, df = vals_l - 1) * sd_boot
  CI_boot_3 <- quantile(boot_sample_means, probs = c(0.025, 0.975), names = FALSE)  #bootstrap 95% CI
  return(list(CI_boot_1 = CI_boot_1, CI_boot_2 = CI_boot_2, CI_boot_3 = CI_boot_3))
}

# matrix to store final results
N_seq_l <- length(N_seq)
out     <- matrix(NA, nrow = N_seq_l, ncol = 5)

set.seed(1)

t_all_1 <- Sys.time()

# iterate over experiment param
for (i in 1:N_seq_l){ # i <- 2; j <- 1
  
  # pull N_tmp (number of observations) for current loop iter
  N_tmp <- N_seq[i]
  ti_1 <- Sys.time()

  # matrix to store results for current N 
  rm(out_tmp)
  out_tmp <-  matrix(NA, nrow = R, ncol = 5)
  
  for (j in 1:R){
    # simulate values N_tmp values from from N(0,1)
    vals      <- rnorm(n = N_tmp, mean = 0, sd = 1)
    vals_mean <- mean(vals)
    vals_sd   <- sd(vals)
    # define CI 
    CI_1 <- vals_mean + c(-1, 1) * qnorm(1 - 0.025) * vals_sd / sqrt(N_tmp)
    CI_2 <- vals_mean + c(-1, 1) * qt(1 - 0.025, df = N_tmp - 1) * vals_sd / sqrt(N_tmp)
    boot_CI_set <- get_boot_CI(vals)
    CI_boot_1 <- boot_CI_set$CI_boot_1
    CI_boot_2 <- boot_CI_set$CI_boot_2
    CI_boot_3 <- boot_CI_set$CI_boot_3
    # get and save CI coverage
    out_tmp[j, ] <- c(
      CI_1[1] * CI_1[2] < 0,
      CI_2[1] * CI_2[2] < 0,
      CI_boot_1[1] * CI_boot_1[2] < 0,
      CI_boot_2[1] * CI_boot_2[2] < 0,
      CI_boot_3[1] * CI_boot_3[2] < 0
      ) * 1
    # remove used elements
    rm(vals, vals_mean, vals_sd, CI_1, CI_2, boot_CI_set, CI_boot_1, CI_boot_2, CI_boot_3)
  }
  # aggregate and save results for current N 
  out[i, ] <- apply(out_tmp, 2, mean)
  
  # prepare final results df 
  out_S <- as.data.frame(out)
  names(out_S) <- c(paste0("ci_approx_", 1:2), paste0("ci_boot_", 1:3))
  out_S$R <- R
  out_S$N <- N_seq
  # save to file
  fwrite(as.data.table(out_S), out_fpath)
  
  # time info
  ti_2 <- Sys.time()
  ti_diff <- round(as.numeric(ti_2 - ti_1, units = "secs"), 3)
  message(paste0("N_tmp: ", N_tmp, " completed in time [s]: ", ti_diff))
}

t_all_2 <- Sys.time()
t_all_diff <- round(as.numeric(t_all_2 - t_all_1, units = "secs"), 3)
message(paste0("All completed in time [s]: ", t_all_diff))
# All completed in time [s]: 2822.672

# prepare final results df 
out <- as.data.frame(out)
names(out_S) <- c(paste0("ci_approx_", 1:2), paste0("ci_boot_", 1:3))
out$R <- R
out$N <- N_seq
# save to file
fwrite(as.data.table(out), out_fpath)


