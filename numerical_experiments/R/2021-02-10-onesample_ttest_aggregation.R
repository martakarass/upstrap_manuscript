
#' This script generate
#' - data results
#' - figures
#' showing how sampling distribution of sample mean and sample standard deviation
#' affect the 
#' 

rm(list = ls())

library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# dir to save results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-02-10-onesample_ttest_aggegating_comparison")

# simulation parameters
N0     <- 50
N1_min <- N0
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu <- 0.3
sigma2_v <- 1
sigma_v  <- 1
rep_R <- 10000


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# get theoretical results for sampling distribution of s sample estimator of sd 
# E(s) = \sigma * c_{4}(n)
# sd(s) = \sigma * \sqrt{1-c_{4}^{2}(n)}
s_mean    <- get_mean_s(n = N0, sigma = sigma_v)
s_sd      <- get_sd_s(n = N0, sigma = sigma_v)
xbar_mean <- mu
xbar_sd   <- sigma_v/sqrt(N0)

# simulate draws from sampling distribution of sample mean 
set.seed(1)
xbar_obs_vec   <- rnorm(rep_R, mean = mu, sd = sigma_v/sqrt(N0))
xbar_obs_vec_l <- length(xbar_obs_vec)

# simulate draws from sampling distribution of sample standard deviation (s) 
set.seed(1)
s_obs_vec   <- rnorm(rep_R, mean = s_mean, sd = s_sd)
s_obs_vec_l <- length(s_obs_vec)

# factor levels, labels 
aggstat_name_level <- c("power_est_median", "power_est_mean")
aggstat_name_label <- c("median", "mean")


# ------------------------------------------------------------------------------
# (1) variate the sample mean estimator 

t1 <- Sys.time()

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_obs_vec[i], sd = s_mean, 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_1 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_1_agg <- 
  plt_df_1 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))


# ------------------------------------------------------------------------------
# (2) variate the sample standard deviation estimator 

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_mean, sd = s_obs_vec[i], 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_2 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_2_agg <- 
  plt_df_2 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))


# ------------------------------------------------------------------------------
# (3) variate the sample mean, (b) sample standard deviation 

mat_out_powerttest <- matrix(NA, nrow = xbar_obs_vec_l, ncol = N1_grid_l)
for (i in 1 : xbar_obs_vec_l){
  mat_out_powerttest[i, ] <- sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = xbar_obs_vec[i], sd = s_obs_vec[i], 
                             sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
}

plt_df_3 <- mat_out_powerttest %>% 
  as.data.frame() %>%
  rename_all(~as.character(N1_grid)) %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols = -idx)%>%
  rename(N1 = name, power_est = value) %>%
  mutate(N1 = as.numeric(N1))
plt_df_3_agg <- 
  plt_df_3 %>%
  group_by(N1) %>%
  summarise(power_est_mean = mean(power_est),
            power_est_median = median(power_est)) %>%
  pivot_longer(cols = c(power_est_mean, power_est_median)) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))

t2 <- Sys.time()
t2 - t1
# 5 min

# sample the idx we store for the purpose of presentation
idx_sub <- seq(1, to = rep_R, length.out = 100)
plt_df_1_sub <- plt_df_1 %>% filter(idx %in% idx_sub)
plt_df_2_sub <- plt_df_2 %>% filter(idx %in% idx_sub)
plt_df_3_sub <- plt_df_3 %>% filter(idx %in% idx_sub)
saveRDS(plt_df_1_sub, file = paste0(results_dir, "/plt_df_1_sub.rds"))
saveRDS(plt_df_2_sub, file = paste0(results_dir, "/plt_df_2_sub.rds"))
saveRDS(plt_df_3_sub, file = paste0(results_dir, "/plt_df_3_sub.rds"))
saveRDS(plt_df_1_agg, file = paste0(results_dir, "/plt_df_1_agg.rds"))
saveRDS(plt_df_2_agg, file = paste0(results_dir, "/plt_df_2_agg.rds"))
saveRDS(plt_df_3_agg, file = paste0(results_dir, "/plt_df_3_agg.rds"))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# plot

plt_df_1_sub <- readRDS(paste0(results_dir, "/plt_df_1_sub.rds"))
plt_df_2_sub <- readRDS(paste0(results_dir, "/plt_df_2_sub.rds"))
plt_df_3_sub <- readRDS(paste0(results_dir, "/plt_df_3_sub.rds"))
plt_df_1_agg <- readRDS(paste0(results_dir, "/plt_df_1_agg.rds"))
plt_df_2_agg <- readRDS(paste0(results_dir, "/plt_df_2_agg.rds"))
plt_df_3_agg <- readRDS(paste0(results_dir, "/plt_df_3_agg.rds"))

# get baseline power 
plt_df_baseline <- data.frame(
  N1 = N1_grid,
  value = sapply(N1_grid, function(n_tmp){
    test_out <- power.t.test(n = n_tmp, delta = 0.3, sd = 1, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
    return(test_out$power)
  })
)

# plot
gg_base_size <- 12

plt1 <- 
  ggplot(plt_df_1_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_1_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $Z_i^{\\bar{X}_N}$, sd = $\\sigma$)'),
       y = TeX('$power.t.test(delta = \\bar{x})$'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  # theme(legend.position = c(0.8, 0.2)) + 
  theme(legend.position = "none",
        legend.title = element_text(size = gg_base_size),
        plot.title = element_text(size = gg_base_size)) + 
  scale_y_continuous(limits = c(0, 1))
# plt1

plt2 <- 
  ggplot(plt_df_2_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_2_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $\\mu$, sd = $Z_i^{S_N}$)'),
       y = TeX('$power.t.test(sd = s)$'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  # theme(legend.position = c(0.8, 0.2)) + 
  theme(legend.position = "none",
        legend.title = element_text(size = gg_base_size),
        plot.title = element_text(size = gg_base_size)) + 
  scale_y_continuous(limits = c(0, 1))
# plt2

plt3 <- 
  ggplot(plt_df_3_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_3_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $Z_i^{\\bar{X}_N}$, sd = $Z_i^{S_N}$)'),
       # title = TeX('power.t.test(delta = $(\\bar{X}_N)_i$, sd = $(S_N)_i$)'),
       y = TeX('$power.t.test(delta = \\bar{x},\\;sd = s$)'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  theme(legend.position = c(0.8, 0.3),
        legend.title = element_text(size = gg_base_size),
        plot.title = element_text(size = gg_base_size)) + 
  scale_y_continuous(limits = c(0, 1))
# plt3

# combined plot
plt_list <- list(plt1, plt2, plt3)
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v")
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.3)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# some further experiments showing the power.t.test values 

# theoretical E(), sd() of sampling distributions for Xbar, S random variable
xbar_mean <- mu
xbar_sd   <- sigma_v/sqrt(N0)
s_mean    <- get_mean_s(n = N0, sigma = sigma_v)
s_sd      <- get_sd_s(n = N0, sigma = sigma_v)

# grid values for the two random variables to plot the density plot
xbar_grid <- seq(
  from = qnorm(p = 0.01, mean = xbar_mean, sd = xbar_sd),
  to = qnorm(p = 1-0.01, mean = xbar_mean, sd = xbar_sd),
  length.out = 10000
)
s_grid <- seq(
  from = qnorm(p = 0.01, mean = s_mean, sd = s_sd),
  to = qnorm(p = 1-0.01, mean = s_mean, sd = s_sd),
  length.out = 10000
)

# combine grids into dfs
## xbar 
arg_power = sapply(xbar_grid, function(xbar_tmp){
  test_out <- power.t.test(n = 100, delta = xbar_tmp, sd = 1, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
  return(test_out$power)
})
plt_xbar_df <- data.frame(
  arg_grid    = xbar_grid,
  arg_density = dnorm(xbar_grid, mean = xbar_mean, sd = xbar_sd),
  arg_power   = arg_power
)
## s 
arg_power = sapply(s_grid, function(s_tmp){
  test_out <- power.t.test(n = 100, delta = 0.3, sd = s_tmp, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
  return(test_out$power)
})
plt_s_df <- data.frame(
  arg_grid    = s_grid,
  arg_density = dnorm(s_grid, mean = s_mean, sd = s_sd),
  arg_power   = arg_power
)

# power draws plot
xbar_draws <- rnorm(n = 10000, mean = xbar_mean, sd = xbar_sd)
s_draws    <- rnorm(n = 10000, mean = s_mean, sd = s_sd)
xbar_power_vals <- sapply(xbar_draws, function(xbar_tmp){
  test_out <- power.t.test(n = 100, delta = xbar_tmp, sd = 1, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
  return(test_out$power)
})
s_power_vals <- sapply(s_draws, function(s_tmp){
  test_out <- power.t.test(n = 100, delta = 0.3, sd = s_tmp, sig.level = 0.05, type = "one.sample", alternative = "two.sided")
  return(test_out$power)
})

plt_xbar_df_agg <- data.frame(
  x = c(mean(xbar_power_vals), median(xbar_power_vals)),
  name = c("power_est_mean", "power_est_median")) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))

plt_s_df_agg <- data.frame(
  x = c(mean(s_power_vals), median(s_power_vals)),
  name = c("power_est_mean", "power_est_median")) %>%
  mutate(name = factor(name, levels = aggstat_name_level, labels = aggstat_name_label))





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# plots 

# density 
plt1a <- 
  ggplot(plt_xbar_df, aes(x = arg_grid, y = arg_density)) + 
  geom_line() + 
  labs(
    # y = TeX('$\\overset{\\bar{X}\\;sampling\\;distribution}{density}$'),
    y = TeX('$\\bar{X}$ distribution density'),
    x = TeX('$\\bar{x}$'))
plt1b <- 
  ggplot(plt_s_df, aes(x = arg_grid, y = arg_density)) + 
  geom_line() + 
  labs(
    y = TeX('$S$ distribution density'),
    x = TeX('$s$'))


# power curve
plt2a <- 
  ggplot(plt_xbar_df, aes(x = arg_grid, y = arg_power)) + 
  geom_line() + 
  # scale_color_futurama() + 
  scale_y_continuous(limits = c(0, 1))  + 
  labs(
    # y = TeX('$\\overset{power.t.test}{n=100,\\;delta = \\bar{x},\\;sd = E\\[S\\]}$'),
    # y = TeX('power.t.test(100, delta = $\\bar{x}$, sd = E\\[S\\])'),
    y = TeX('$\\phi(\\bar{x}) = power.t.test(delta = \\bar{x}$)'),
    color = "power.t.test\nvalues\naggregate",
    x = TeX('$\\bar{x}$')) 
  # theme(legend.position = c(0.8, 0.3))
plt2b <- 
  ggplot(plt_s_df, aes(x = arg_grid, y = arg_power)) + 
  geom_line() + 
  scale_color_futurama() + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(
    y = TeX('$\\phi(s) = power.t.test(sd = s$)'),
    color = "power.t.test\nvalues\naggregate",
    x = TeX('$s$'))

# power density
plt3a <- 
  ggplot(data.frame(x = xbar_power_vals), aes(x = x)) + 
  geom_density() + 
  labs(
    # y = TeX('$\\overset{\\bar{X}\\;sampling\\;distribution}{density}$'),
    y = TeX('$\\phi(\\bar{X})$ distribution density'),
    x = TeX('$\\phi(\\bar{x})$'),
    color =  TeX('$\\phi(\\bar{X})$: ')) + 
  scale_x_continuous(limits = c(0, 1)) +   
  geom_vline(data = plt_xbar_df_agg, aes(xintercept = x, color = name)) + 
  scale_color_futurama() + 
  theme(legend.position = c(0.2, 0.8))
plt3b <- 
  ggplot(data.frame(x = s_power_vals), aes(x = x)) + 
  geom_density() + 
  labs(
    # y = TeX('$\\overset{\\bar{X}\\;sampling\\;distribution}{density}$'),
    y = TeX('$\\phi(S)$ distribution density'),
    x = TeX('$\\phi(s)$'),
    color =  TeX('$\\phi(S)$: ')) + 
  scale_x_continuous(limits = c(0, 1)) +   
  geom_vline(data = plt_s_df_agg, aes(xintercept = x, color = name)) + 
  scale_color_futurama() + 
  theme(legend.position = c(0.2, 0.8))


# combined plot
plt_list <- list(plt1a, plt1b, 
                 plt2a, plt2b, 
                 plt3a, plt3b)
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = FALSE)
plt


plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_2.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 6)


