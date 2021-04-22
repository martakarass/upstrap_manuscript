

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)
source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))


# ------------------------------------------------------------------------------
# PLOT 1 

# read simulated example data
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
out_df <- readRDS(paste0(results_dir, "/out_df.rds"))

# aggregate data 
out_df_agg <- out_df %>%
  group_by(M_sample_size) %>%
  summarise(
    powerttest_mean = mean(powerttest),
    powerttest_median = median(powerttest),
    upstrap_mean = mean(upstrap),
    upstrap_median = median(upstrap),
    ttest_mean = mean(ttest)
  ) %>%
  pivot_longer(cols = -M_sample_size) %>%
  separate(name, into = c("est_approach", "results_agg_method"))
head(out_df_agg)

# define gold standard power data 
n_arg <- sort(unique(out_df$M_sample_size))
goldstandard_df <- data.frame(
  M_sample_size = n_arg,
  value = power.t.test(n = n_arg, delta = 0.3, sd = 1, type = "one.sample")$power
)

# pivot pre-aggregated data to generate a plot
out_df_L <- 
  out_df %>%
  pivot_longer(cols = -c(rep_idx, M_sample_size), names_to = "est_approach") 
head(out_df_L)

est_approach_levels <- c("powerttest", "upstrap", "ttest")
est_approach_labels <- c("(a) power.t.test", "(b) upstrap", "(c) t.test")

plt_list <- list()
for (i in 1 : length(est_approach_levels)){ # i <- 2
  est_approach_i  <- est_approach_levels[i]
  out_df_agg_i    <- out_df_agg %>% filter(est_approach == est_approach_i)
  out_df_L_i      <- out_df_L %>% filter(est_approach == est_approach_i)
  plt <- ggplot()
  if (i %in% c(1,2)){
    plt <- plt +
      geom_line(data = out_df_L_i %>% filter(rep_idx <= 100),
                aes(x = M_sample_size, y = value, group = rep_idx),
                color = "black", alpha = 0.1, size = 0.3) 
  }
  plt <- plt + 
    geom_line(data = out_df_agg_i, 
              aes(x = M_sample_size, y = value, color = results_agg_method, group = results_agg_method),
              size = 1, inherit.aes = FALSE) + 
    geom_line(data = goldstandard_df, aes(x = M_sample_size, y = value, group = 1), 
              color = "black", inherit.aes = FALSE) + 
    theme(
      legend.title=element_text(size = legend_font_size_numexp),
      plot.title = element_text(size = title_font_size_numexp),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = facetrid_font_size_numexp)
    ) + 
    labs(x = TeX("Sample size $M_k$"),
         y = "",
         color = "",
         title = est_approach_labels[i]) + 
    scale_y_continuous(limits = c(0,1))  + 
    theme(legend.position = c(0.7, 0.3)) 
  # append plot to plots list
  plt_list[[length(plt_list) + 1]] <- plt 
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = TRUE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_PLOT1.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.5)


# ------------------------------------------------------------------------------
# PLOT 2

rm(out_df, out_df_L, out_df_agg, goldstandard_df)

# read simulated example data
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
out_df <- readRDS(paste0(results_dir, "/out_df_PLOT2.rds"))
head(out_df)

# aggregate data 
out_df_agg <- out_df %>%
  group_by(M_sample_size, example) %>%
  summarise(
    powerttest_mean = mean(powerttest),
    powerttest_median = median(powerttest)
  ) %>%
  pivot_longer(cols = -c(M_sample_size, example)) %>%
  separate(name, into = c("est_approach", "results_agg_method"))
head(out_df_agg)

# define gold standard power data 
n_arg <- sort(unique(out_df$M_sample_size))
goldstandard_df <- data.frame(
  M_sample_size = n_arg,
  value = power.t.test(n = n_arg, delta = 0.3, sd = 1, type = "one.sample")$power
)


example_levels <- c("ex_a", "ex_b", "ex_c")
example_labels <- c('$power.t.test(delta = \\bar{x},\\;sd = 1$)', 
                    '$power.t.test(delta = 0.3,\\;sd = s$)', 
                    '$power.t.test(delta = \\bar{x},\\;sd = s$)')

plt_list <- list()
for (i in 1 : length(example_levels)){ # i <- 1
  example_i     <- example_levels[i]
  out_df_i      <- out_df %>% filter(example == example_i)
  out_df_agg_i  <- out_df_agg %>% filter(example == example_i)
  plt <- 
    ggplot(out_df_i %>% filter(rep_idx <= 100), 
           aes(x = M_sample_size, y = powerttest, group = rep_idx)) + 
    geom_line(alpha = 0.1, size = 0.3) + 
    geom_line(data = out_df_agg_i, 
              aes(x = M_sample_size, y = value, color = results_agg_method, group = results_agg_method),
              size = 1) + 
    geom_line(data = goldstandard_df, aes(x = M_sample_size, y = value, group = 1), 
              color = "black") + 
    theme(
      legend.title=element_text(size = legend_font_size_numexp),
      plot.title = element_text(size = title_font_size_numexp),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = facetrid_font_size_numexp)
    ) + 
    labs(x = TeX("Sample size $M_k$"),
         y = "",
         color = "",
         title = TeX(example_labels[i])) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(legend.position = c(0.7, 0.3)) 
  # append plot to plots list
  plt_list[[length(plt_list) + 1]] <- plt 
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = TRUE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_PLOT2.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.5)



# ------------------------------------------------------------------------------
# PLOT 3

rm(out_df, out_df_L, out_df_agg, goldstandard_df)

# read simulated example data
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-04-05-onesample_ttest_aggegating_comparison")
out_df <- readRDS(paste0(results_dir, "/out_df_PLOT3.rds"))
head(out_df)

plt1 <- 
  out_df %>% 
  filter(arg_name == "xbar") %>%
  ggplot(aes(x = arg_val)) + 
  geom_density()
plt1


# simulation parameters
N   <- 50
R   <- 100000
mu  <- 0.3
sigma2 <- 1
M_k <- 100

xbar_obs_vec <- seq(-0.1, 0.7, length.out = 1000)
s_obs_vec <- seq(0.75, 1.25, length.out = 1000)



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

# plot final: combined plot
plt_list <- list(plt1a, plt1b, 
                 plt2a, plt2b)
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = FALSE)
plt


plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_PLOT3.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 6)









