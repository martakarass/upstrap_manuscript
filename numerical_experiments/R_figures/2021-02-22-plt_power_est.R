
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

my_pal <- pal_futurama()(12)
show_col(my_pal)
show_col(my_pal[c(4,9)])


# PLOT 1a: one-sample t-test power  ---------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-onesample_ttest_agg/")
onesample_power_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_power_est_bootCI.rds"))
onesample_samplesize_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_samplesize_est_bootCI.rds"))
onesample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret.rds"))

N0_levels <- sort(unique(onesample_power_df$N0))
plt_df <- 
  onesample_power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_1a <- 
  ggplot(
    plt_df, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_line() + 
  geom_line(data = onesample_power_theoret_df, aes(x = N1, y = power_est, group =1), inherit.aes = FALSE, alpha = 0.3, size = 1) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(plt_df_N_pointer$N0, 100, 150, 200)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Problem: one-sample t-test") + 
  scale_color_manual(values = my_pal[c(4,9)]) 
plt_1a


# PLOT 1b: two-sample t-test power  ---------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-twosample_ttest_agg")
twosample_power_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_power_est_bootCI.rds"))
twosample_samplesize_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_samplesize_est_bootCI.rds"))
twosample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret"))

N0_levels <- sort(unique(twosample_power_df$N0))
plt_df <- 
  twosample_power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_1b <- 
  ggplot(
    plt_df, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_line() + 
  geom_line(data = twosample_power_theoret_df, aes(x = N1, y = power_est, group =1), inherit.aes = FALSE, alpha = 0.3, size = 1) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(plt_df_N_pointer$N0, 200, 250, 300)) + 
  theme(
    legend.position = c(0.8, 0.1),
    legend.title=element_text(size = 12),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Problem: two-sample t-test",
       color = "Upstrap type:") + 
  scale_color_manual(values = my_pal[c(4,9)]) 
plt_1b


# combine the two plots
plt_list = list(plt_1a, plt_1b)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/ttest_power_est.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 8)



# PLOT 1b: two-sample t-test power  ---------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-twosample_ttest_agg")
twosample_power_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_power_est_bootCI.rds"))
twosample_samplesize_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_samplesize_est_bootCI.rds"))
twosample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret"))

N0_levels <- sort(unique(twosample_power_df$N0))
plt_df <- 
  twosample_power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_1b <- 
  ggplot(
    plt_df, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_line() + 
  geom_line(data = twosample_power_theoret_df, aes(x = N1, y = power_est, group =1), inherit.aes = FALSE, alpha = 0.3, size = 1) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(plt_df_N_pointer$N0, 200, 250, 300)) + 
  theme(
    legend.position = c(0.8, 0.1),
    legend.title=element_text(size = 12),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Problem: two-sample t-test",
       color = "Upstrap type:") + 
  scale_color_manual(values = my_pal[c(4,9)]) 
plt_1b


# combine the two plots
plt_list = list(plt_1a, plt_1b)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/ttest_power_est.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 8)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 2a: LMM

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-lmm_trt_agg")
power_df <- readRDS(paste0(fdir_tmp, "/lmm_trt_power_est_bootCI.rds"))

N0_levels <- sort(unique(power_df$N0))
plt_df_boot <- 
  power_df %>% 
  filter(method_name == "out_boot_GEE") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_gold <- 
  power_df %>% 
  filter(method_name == "out_gold_GEE") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_boot_N_pointer <- plt_df_boot %>% filter(N0 == N1)

plt_2a <- 
  ggplot(
    plt_df_boot, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line() + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group =1), inherit.aes = FALSE, alpha = 0.4, size = 0.7) + 
  geom_ribbon(data = plt_df_gold, aes(ymin = est_bootci_mean_lwr, ymax = est_bootci_mean_upr, fill = type), 
              alpha = 0.2, color = "grey50", fill = NA, size = 0.5, linetype = 3) + 
  
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(0, plt_df_boot_N_pointer$N0, 50, 100, 150)) +
  theme(
    legend.position = "none",
    # legend.position = c(0.8, 0.1),
    # legend.title=element_text(size = 12),
    # legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Problem: multilevel linear regression (GEE est. method)",
       color = "Upstrap type:") + 
  scale_color_manual(values = my_pal[c(4,9)]) + 
  scale_fill_manual(values = my_pal[c(4,9)])
plt_2a



# ------------------------------------------------------------------------------
# PLOT 2b: GLMM

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-glmm_trt_agg")
power_df <- readRDS(paste0(fdir_tmp, "/glmm_trt_power_est_bootCI.rds"))

N0_levels <- sort(unique(power_df$N0))
plt_df_boot <- 
  power_df %>% 
  filter(method_name == "out_boot_GEE") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_gold <- 
  power_df %>% 
  filter(method_name == "out_gold_GEE") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_boot_N_pointer <- plt_df_boot %>% filter(N0 == N1)

plt_2b <- 
  ggplot(
    plt_df_boot, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line() + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group =1), inherit.aes = FALSE, alpha = 0.4, size = 0.7) + 
  geom_ribbon(data = plt_df_gold, aes(ymin = est_bootci_mean_lwr, ymax = est_bootci_mean_upr, fill = type), 
              alpha = 0.2, color = "grey50", fill = NA, size = 0.5, linetype = 3) + 
  
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(plt_df_boot_N_pointer$N0, 150, 200, 250)) +
  theme(
    legend.position = c(0.8, 0.25),
    legend.title=element_text(size = 12),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12)
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Problem: multilevel generalized linear regression (GEE est. method)",
       color = "Upstrap type:", fill = "Upstrap type:") + 
  scale_color_manual(values = my_pal[c(4,9)]) + 
  scale_fill_manual(values = my_pal[c(4,9)])
plt_2b


# combine the two plots
plt_list = list(plt_2a, plt_2b)
plt <- plot_grid(plotlist = plt_list, ncol = 1, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/lmm_glmm_power_est.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 7, base_height = 10)
