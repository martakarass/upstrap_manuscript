
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

show_col(my_pal)
my_pal_numexp2 <- my_pal[c(7,8)]
show_col(my_pal_numexp)
show_col(my_pal_numexp2)


# PLOT 1 -----------------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-05-25-onesample_ttest_compare_upstrap_powerttest_agg/")
onesample_power_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_compare_upstrap_powerttest_bootCI.rds"))

N0_levels <- sort(unique(onesample_power_df$N0))
plt_df <- 
  onesample_power_df %>% 
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)

range(plt_df$N1)
plt_x_axis_max <- 150

# trick to handle Error: Aesthetics can not vary with a ribbon, so we duplicate
# the bridge observation so as it has both levels
plt_df_dupl <- plt_df %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df)

ylim_max <- 0.015

plt_1 <- 
  ggplot(
    plt_df_dupl, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), 
              alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # geom_line(data = onesample_power_theoret_df, aes(x = N1, y = power_est, group =1), 
  #           inherit.aes = FALSE, size = 0.4) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) +
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, alpha = 1, size = 0.2, color = "black") +
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(-1, 1) * ylim_max, breaks = seq(-ylim_max, ylim_max, by = 0.005)) +
  # scale_x_continuous(breaks = c(0, plt_df_N_pointer$N0, 100, 150), limits = c(0, 150)) + 
  scale_x_continuous(breaks = seq(0, 150, by = 50), limits = c(0, 150)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power difference", 
       title = "Problem 1: linear regression with intercept only\n(Est. method: one-sample t-test)") + 
  scale_color_manual(values = my_pal_numexp2) + 
  scale_fill_manual(values = my_pal_numexp2) 
plt_1



# PLOT 2 -----------------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-05-25-twosample_ttest_compare_upstrap_powerttest_agg/")
onesample_power_df <- readRDS(paste0(fdir_tmp, "/twosample_ttest_compare_upstrap_powerttest_bootCI.rds"))

N0_levels <- sort(unique(onesample_power_df$N0))
plt_df <- 
  onesample_power_df %>% 
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)

range(plt_df$N1)
plt_x_axis_max <- 300

plt_df_dupl <- plt_df %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df)

# ylim <- c(-1, 1) * max(abs(c(plt_df_dupl$power_est_aggmedian_lwr, plt_df_dupl$power_est_aggmedian_upr))) 

plt_2 <- 
  ggplot(
    plt_df_dupl, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), 
              alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # geom_line(data = onesample_power_theoret_df, aes(x = N1, y = power_est, group =1), 
  #           inherit.aes = FALSE, size = 0.4) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) +
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, alpha = 1, size = 0.2, color = "black") +
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(-1, 1) * ylim_max, breaks = seq(-ylim_max, ylim_max, by = 0.005)) +
  # scale_x_continuous(breaks = c(0, plt_df_N_pointer$N0, 100, 150), limits = c(0, 150)) + 
  scale_x_continuous(breaks = seq(0, 300, by = 50), limits = c(0, 300)) + 
  theme(
    legend.position = c(0.8, 0.08),
    legend.title=element_text(size = legend_font_size_numexp),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power difference", 
       title = "Problem 2: simple linear regression\n(Est. method: two-sample t-test)",
       color = "Upstrap type:",
       fill = "Upstrap type:"
       ) + 
  scale_color_manual(values = my_pal_numexp2) + 
  scale_fill_manual(values = my_pal_numexp2) 
plt_2


# combine the two plots
plt_list = list(plt_1, plt_2)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/ttest_compare_upstrap_and_powerttest.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 12)





