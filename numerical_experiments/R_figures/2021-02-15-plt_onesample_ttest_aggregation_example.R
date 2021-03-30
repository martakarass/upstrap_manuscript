
#' This script generate
#' - figures
#' showing how sampling distribution of sample mean and sample standard deviation
#' affect the (a) mean, (b) median of power across multiple experiment repetitions. 

rm(list = ls())

library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# dir to results 
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-02-10-onesample_ttest_aggegating_comparison")

# # simulation parameters
N0     <- 50
N1_min <- N0
N1_max <- 200
N1_grid <- N1_min : N1_max
# N1_grid_l <- length(N1_grid)
# 
# # data generating model
mu <- 0.3
sigma2_v <- 1
sigma_v  <- 1
# rep_R <- 10000

# factor levels, labels 
aggstat_name_level <- c("power_est_median", "power_est_mean")
aggstat_name_label <- c("median", "mean")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PART 1: power curves and their aggregates: (a) mean, (b) median
 
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

plt1 <- 
  ggplot(plt_df_1_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_1_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size_numexp) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $Z_i^{\\bar{X}_N}$, sd = $\\sigma$)'),
       y = TeX('$power.t.test(delta = \\bar{x})$'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  # theme(legend.position = c(0.8, 0.2)) + 
  theme(legend.position = "none",
        legend.title = element_text(size = legend_font_size_numexp),
        # plot.title = element_text(size = gg_base_size_numexp)
        ) + 
  scale_y_continuous(limits = c(0, 1))
# plt1

plt2 <- 
  ggplot(plt_df_2_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_2_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size_numexp) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $\\mu$, sd = $Z_i^{S_N}$)'),
       y = TeX('$power.t.test(sd = s)$'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  # theme(legend.position = c(0.8, 0.2)) + 
  theme(legend.position = "none",
        legend.title = element_text(size = legend_font_size_numexp),
        # plot.title = element_text(size = gg_base_size_numexp)
  ) + 
  scale_y_continuous(limits = c(0, 1))
# plt2

plt3 <- 
  ggplot(plt_df_3_sub,
         aes(x = N1, y = power_est, group = idx)) + 
  geom_line(alpha = 0.2, size = 0.1) + 
  geom_line(data = plt_df_baseline, aes(x = N1, y = value, group = 1), inherit.aes = FALSE, size = 1.3) + 
  geom_line(data = plt_df_3_agg, aes(x = N1, y = value, color = name), inherit.aes = FALSE) + 
  theme_classic(base_size = gg_base_size_numexp) + 
  labs(x = "Sample size", 
       # y = "Power",
       # title = TeX('power.t.test(delta = $Z_i^{\\bar{X}_N}$, sd = $Z_i^{S_N}$)'),
       # title = TeX('power.t.test(delta = $(\\bar{X}_N)_i$, sd = $(S_N)_i$)'),
       y = TeX('$power.t.test(delta = \\bar{x},\\;sd = s$)'),
       color = "point-wise\naggregate") + 
  scale_color_futurama() + 
  theme(legend.position = c(0.8, 0.3),
        legend.title = element_text(size = legend_font_size_numexp),
        # plot.title = element_text(size = gg_base_size_numexp)
  ) + 
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
# PART 2: dive into function: power(xbar_i), power(s_i), power(xbar_i, s_i)

# get theoretical E(), sd() of sampling distributions for Xbar, S random variables
xbar_mean <- mu
xbar_sd   <- sigma_v/sqrt(N0)
s_mean    <- get_mean_s(n = N0, sigma = sigma_v)
s_sd      <- get_sd_s(n = N0, sigma = sigma_v)

# grid values for Xbar, S random variables to plot the density plot
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

# get "draws" from power 
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
# generate the plots 

# Plot 1: density Xbar, density S
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


# Plot 2
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

# Plot 3: power density
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


# plot final: combined plot
plt_list <- list(plt1a, plt1b, 
                 plt2a, plt2b, 
                 plt3a, plt3b)
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = FALSE)
plt


plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_2.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 6)


