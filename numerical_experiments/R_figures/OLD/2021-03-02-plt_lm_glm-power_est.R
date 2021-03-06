
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

my_pal <- pal_futurama()(12)
my_pal_sub <- my_pal[c(4,9)]

show_col(my_pal)
show_col(my_pal_sub)


# PLOT 1 -----------------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-ttest_agg/")
power_df <- readRDS(paste0(fdir_tmp, "/ttest_power_est_bootCI.rds"))
power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret.rds"))

N0_levels <- sort(unique(power_df$N0))
plt_df <- 
  power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_x_axis_max <- plt_df %>% filter(N0 == max(N0)) %>% filter(power_est_aggmedian >= 0.95) %>% pull(N1) %>% min()
plt_x_axis_max

plt_1 <- 
  ggplot(
    plt_df, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line() + 
  geom_line(data = power_theoret_df, aes(x = N1, y = power_est, group =1), inherit.aes = FALSE, alpha = 0.5) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(0, plt_df_N_pointer$N0, 100, 150), limits = c(0, 150)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank()
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "One-sample t-test") + 
  scale_color_manual(values = my_pal_sub) + 
  scale_fill_manual(values = my_pal_sub) 
plt_1



# PLOT 2 -----------------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-twosample_ttest_agg")
twosample_power_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_power_est_bootCI.rds"))
twosample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret"))

N0_levels <- sort(unique(twosample_power_df$N0))
plt_df <- 
  twosample_power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_x_axis_max <- plt_df %>% filter(N0 == max(N0)) %>% filter(power_est_aggmedian >= 0.95) %>% pull(N1) %>% min()
plt_x_axis_max

plt_2 <- 
  ggplot(
    plt_df, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line() + 
  geom_line(data = twosample_power_theoret_df, aes(x = N1, y = power_est, group =1), inherit.aes = FALSE, alpha = 0.5) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = c(plt_df_N_pointer$N0, 200, 250, 300)) + 
  theme(
    legend.position = c(0.8, 0.1),
    legend.title=element_text(size = 10),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = 12),
    panel.grid.major.x = element_blank()
  ) + 
  labs(x = "Sample size M", y = "Estimated power", title = "Two-sample t-test",
       color = "Upstrap type:",
       fill = "Upstrap type:"
       ) + 
  scale_color_manual(values = my_pal_sub) + 
  scale_fill_manual(values = my_pal_sub) 
plt_2


# combine the two plots
plt_list = list(plt_1, plt_2)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

# plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/2021-03-01-ttest_power_est.png")
plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/ttest_power_est.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 6)





