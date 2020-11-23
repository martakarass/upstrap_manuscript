
rm(list = ls())

# libraries
library(data.table)
library(tidyverse)
# devtools::install_github("ccrainic/upstrap")
# devtools::install_github("martakarass/upstrap")
# library(upstrap)

# parameters definition
project_dir <- "/Users/marta/Dropbox/_PROJECTS/upstrap_manuscript"
mu        <- 0.1
R_boot    <- 1000
power_val <- 0.8


# get gold standard 
n_grid <- 200:1200
power_out <- sapply(n_grid, function(n){
  (power.t.test(n = n, delta = mu, sd = 1, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
})
  
plt_df <- data.frame(n_grid, power_out)
n_grid_whichmin <- plt_df %>% filter(power_out >= 0.8) %>% 
  filter(n_grid == min(n_grid)) %>% pull(n_grid) %>% unlist()
plt <- 
  ggplot(plt_df, aes(x = n_grid, y = power_out)) + 
  geom_hline(yintercept = 0.8, color = "blue") + 
  geom_vline(xintercept = n_grid_whichmin, linetype = 3, color = "red") + 
  geom_line(color = "red", size = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = seq(100, 1200, by = 100)) + 
  labs(x = "Sample size n", y = "power to reject H0") 
plt


# ------------------------------------------------------------------------------
# estimate power with power.t.test starting with some small sample 

rep_n <- 10000
# n0 <- 500
n0 <- 200
n0_n_grid <- n0:1200

out_mat <- matrix(NA, ncol = length(n0_n_grid), nrow = rep_n)

set.seed(1)
for (i in 1:rep_n){
  if (i %% 10 == 0) print(i)
  sample_i <- rnorm(n = n0, mean = mu, sd = 1)
  sample_i_mean <- mean(sample_i)
  sample_i_sd   <- sd(sample_i)
  out_mat[i, ] <- sapply(n0_n_grid, function(n){
    (power.t.test(n = n, delta = sample_i_mean, sd = sample_i_sd, sig.level = 0.05, type = "one.sample", alternative = "two.sided"))$power
  })
}

out_df_agg <- data.frame(
  power_out_mean = apply(out_mat, 2, mean),
  power_out_sd = apply(out_mat, 2, sd),
  power_out_median = apply(out_mat, 2, median),
  n_grid = n0_n_grid,
  n0 = n0
)

n_grid_whichmin <- plt_df %>% filter(power_out >= 0.8) %>% 
  filter(n_grid == min(n_grid)) %>% pull(n_grid) %>% unlist()
plt <- 
  ggplot(plt_df, aes(x = n_grid, y = power_out)) + 
  geom_hline(yintercept = 0.8, color = "blue") + 
  geom_vline(xintercept = n_grid_whichmin, linetype = 3, color = "red") + 
  geom_line(color = "red", size = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  geom_line(data = out_df_agg, aes(x = n_grid, y = power_out_mean), color = "black") + 
  geom_line(data = out_df_agg, aes(x = n_grid, y = power_out_median), color = "brown") + 
  scale_x_continuous(breaks = seq(100, 1200, by = 100)) + 
  labs(x = "Sample size n", y = "power to reject H0")  + 
  theme_gray(base_size = 16)
plt

# out_fpath <-  paste0(project_dir, "/numerical_experiments/results/2020-11-18-ttest_normal_mean/n0_500.rds")
out_fpath <-  paste0(project_dir, "/numerical_experiments/results/2020-11-18-ttest_normal_mean/n0_200.rds")
saveRDS(object = out_df_agg, file = out_fpath)


# ------------------------------------------------------------------------------
# plot the distributions of vals 

nval <- n_grid_whichmin
colidx <- which(n0_n_grid == nval)
plt_df_2 <- data.frame(x = out_mat[, colidx])
x_mean   <- mean(plt_df_2$x)
x_median <- median(plt_df_2$x)

ggplot(plt_df_2, aes(x = x)) + 
  geom_histogram(binwidth = 0.01, fill = "grey60", color = "black") + 
  labs(x = paste0("Power for n=", nval, "\nest. from R=10,000 samples of n0=200"), y = "")  + 
  theme_gray(base_size = 22) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  geom_vline(xintercept = x_median, color = "brown", size = 2) + 
  geom_vline(xintercept = x_mean, color = "black", size = 2)


nval <- 200
colidx <- which(n0_n_grid == nval)
plt_df_2 <- data.frame(x = out_mat[, colidx])
x_mean   <- mean(plt_df_2$x)
x_median <- median(plt_df_2$x)

ggplot(plt_df_2, aes(x = x)) + 
  geom_histogram(binwidth = 0.01, fill = "grey60", color = "black") + 
  labs(x = paste0("Power for n=", nval, "\nest. from R=10,000 samples of n0=200"), y = "")  + 
  theme_gray(base_size = 22) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  geom_vline(xintercept = x_median, color = "brown", size = 2) + 
  geom_vline(xintercept = x_mean, color = "black", size = 2)

nval <- 1200
colidx <- which(n0_n_grid == nval)
plt_df_2 <- data.frame(x = out_mat[, colidx])
x_mean   <- mean(plt_df_2$x)
x_median <- median(plt_df_2$x)

ggplot(plt_df_2, aes(x = x)) + 
  geom_histogram(binwidth = 0.01, fill = "grey60", color = "black") + 
  labs(x = paste0("Power for n=", nval, "\nest. from R=10,000 samples of n0=200"), y = "")  + 
  theme_gray(base_size = 22) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
  geom_vline(xintercept = x_median, color = "brown", size = 2) + 
  geom_vline(xintercept = x_mean, color = "black", size = 2)








