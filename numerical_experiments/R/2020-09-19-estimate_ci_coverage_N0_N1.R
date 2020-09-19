
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' @description 
#' Script to estimate CI coverage of parameter true value in the experiment
#' of estimating N(0,1) sample mean. 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>

## -----------------------------------------------------------------------------

## args
# arg_str <- "N0_100_N1_150_R_50_BBOOT_10000"
arg_str <- as.character(args[1])

## fixed params
#project_dir <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript"
project_dir <- "/users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments"

# derivative args
N0 <- strsplit(arg_str, split = "_")[[1]][2]
N1 <- strsplit(arg_str, split = "_")[[1]][4]
R  <- strsplit(arg_str, split = "_")[[1]][6]
B_boot  <- strsplit(arg_str, split = "_")[[1]][8]
N0 <- as.numeric(N0)
N1 <- as.numeric(N1)
R  <- as.numeric(R)
B_boot  <- as.numeric(B_boot)
out_fname <- paste0(arg_str, ".csv")
out_fpath <- paste0(project_dir,"/numerical_experiments/results_CL/2020-09-19-estimate_ci_coverage_N0_N1_results/", out_fname)
# message with derivative args values
message(paste0(c("N0", "N1", "R", "B_boot"), collapse = ", "))
message(paste0(c(N0, N1, R, B_boot), collapse = ", "))

# libraries
library(data.table)
library(dplyr)


## -----------------------------------------------------------------------------

set.seed(1)

# object to store results from the repetition of the experiment 
mat_out_rrs <- matrix(NA, nrow = R, ncol = 11)

t_all_1 <- Sys.time()

for (rr in 1:R){
  print(rr)
  
  # simulate observed sample
  vals      <- rnorm(n = N0, mean = 0, sd = 1)
  vals_mean <- mean(vals)
  vals_sd   <- sd(vals)
  
  # boot resamples
  boot_samples <- matrix(sample(x = vals, size = (B_boot * N1), replace = TRUE), nrow = B_boot, ncol = N1)
  boot_samples_stat <- apply(boot_samples, 1, mean)
  boot_m  <- mean(boot_samples_stat) 
  boot_sd <- sd(boot_samples_stat) 
  # boot CI
  CI_boot_1 <- boot_m + c(-1, 1) * qnorm(1 - 0.025) * boot_sd 
  CI_boot_2 <- boot_m + c(-1, 1) * qt(1 - 0.025, df = N0 - 1) * boot_sd
  CI_boot_3 <- boot_m + c(-1, 1) * qt(1 - 0.025, df = N1 - 1) * boot_sd
  CI_boot_4 <- quantile(boot_samples_stat, probs = c(0.025, 0.975), names = FALSE)  #bootstrap 95% CI
  # boot CI coverage 
  CI_coverage <- sapply(list(CI_boot_1, CI_boot_2, CI_boot_3, CI_boot_4), function(CI) {(CI[1] * CI[2] < 0) * 1}) 
  CI_length   <- sapply(list(CI_boot_1, CI_boot_2, CI_boot_3, CI_boot_4), function(CI) {diff(CI)}) 
  # append results from the current repetition of the experiment 
  out_rr <- c(
    rr,
    boot_m,
    boot_sd,
    CI_coverage,
    CI_length
  )
  mat_out_rrs[rr, ] <- out_rr
  
  # aggregate and save every other repetition
  if (rr %% 1000 == 0){
    message(paste0("rr: ", rr))
    t_all_2 <- Sys.time()
    t_all_diff <- round(as.numeric(t_all_2 - t_all_1, units = "secs"), 3)
    # aggregate
    df_out_rrs <- as.data.frame(mat_out_rrs, stringsAsFactors = FALSE) 
    names(df_out_rrs) <- c("rr", "boot_m", "boot_sd", paste0("CI_coverage_", 1:4), paste0("CI_length_", 1:4))
    df_out_agg <- cbind(
      df_out_rrs %>% select(starts_with("boot_m"))  %>% summarise_all(mean, na.rm = TRUE)  %>% setNames("boot_m_mean"),
      df_out_rrs %>% select(starts_with("boot_m"))  %>% summarise_all(sd, na.rm = TRUE)    %>% setNames("boot_m_sd"),
      df_out_rrs %>% select(starts_with("boot_sd")) %>% summarise_all(mean, na.rm = TRUE)  %>% setNames("boot_sd_mean"),
      df_out_rrs %>% select(starts_with("boot_sd")) %>% summarise_all(sd, na.rm = TRUE)    %>% setNames("boot_sd_sd"),
      df_out_rrs %>% select(starts_with("CI_coverage_")) %>% summarise_all(mean, na.rm = TRUE) %>% setNames(paste0("CI_coverage_", 1:4, "_mean")),
      df_out_rrs %>% select(starts_with("CI_coverage_")) %>% summarise_all(sd, na.rm = TRUE)   %>% setNames(paste0("CI_coverage_", 1:4, "_sd")),
      df_out_rrs %>% select(starts_with("CI_length_")) %>% summarise_all(mean, na.rm = TRUE) %>% setNames(paste0("CI_length_", 1:4, "_mean")),
      df_out_rrs %>% select(starts_with("CI_length_")) %>% summarise_all(sd, na.rm = TRUE)   %>% setNames(paste0("CI_length_", 1:4, "_sd"))
    )
    df_out_agg$cnt <- sum(!is.na(df_out_rrs$boot_m))
    df_out_agg$cnt <- sum(!is.na(df_out_rrs$boot_m))
    df_out_agg$exec_secs <- t_all_diff
    df_out_agg$N0 <- N0
    df_out_agg$N1 <- N1
    # save
    fwrite(as.data.table(df_out_agg), out_fpath)
  }
  
}

t_all_2 <- Sys.time()
t_all_diff <- round(as.numeric(t_all_2 - t_all_1, units = "secs"), 3)
# aggregate
df_out_rrs <- as.data.frame(mat_out_rrs, stringsAsFactors = FALSE) 
names(df_out_rrs) <- c("rr", "boot_m", "boot_sd", paste0("CI_coverage_", 1:4), paste0("CI_length_", 1:4))
df_out_agg <- cbind(
  df_out_rrs %>% select(starts_with("boot_m"))  %>% summarise_all(mean, na.rm = TRUE)  %>% setNames("boot_m_mean"),
  df_out_rrs %>% select(starts_with("boot_m"))  %>% summarise_all(sd, na.rm = TRUE)    %>% setNames("boot_m_sd"),
  df_out_rrs %>% select(starts_with("boot_sd")) %>% summarise_all(mean, na.rm = TRUE)  %>% setNames("boot_sd_mean"),
  df_out_rrs %>% select(starts_with("boot_sd")) %>% summarise_all(sd, na.rm = TRUE)    %>% setNames("boot_sd_sd"),
  df_out_rrs %>% select(starts_with("CI_coverage_")) %>% summarise_all(mean, na.rm = TRUE) %>% setNames(paste0("CI_coverage_", 1:4, "_mean")),
  df_out_rrs %>% select(starts_with("CI_coverage_")) %>% summarise_all(sd, na.rm = TRUE)   %>% setNames(paste0("CI_coverage_", 1:4, "_sd")),
  df_out_rrs %>% select(starts_with("CI_length_")) %>% summarise_all(mean, na.rm = TRUE) %>% setNames(paste0("CI_length_", 1:4, "_mean")),
  df_out_rrs %>% select(starts_with("CI_length_")) %>% summarise_all(sd, na.rm = TRUE)   %>% setNames(paste0("CI_length_", 1:4, "_sd"))
)
df_out_agg$cnt <- sum(!is.na(df_out_rrs$boot_m))
df_out_agg$cnt <- sum(!is.na(df_out_rrs$boot_m))
df_out_agg$exec_secs <- t_all_diff
df_out_agg$N0 <- N0
df_out_agg$N1 <- N1
# save
fwrite(as.data.table(df_out_agg), out_fpath)


