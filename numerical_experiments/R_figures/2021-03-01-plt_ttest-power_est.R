
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

show_col(my_pal)
show_col(my_pal_numexp)


# PLOT 1 -----------------------------------------------------------------------

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-onesample_ttest_agg/")
onesample_power_df <- readRDS(paste0(fdir_tmp, "/onesample_ttest_power_est_bootCI.rds"))
onesample_power_theoret_df <- readRDS(paste0(fdir_tmp, "/res_theoret.rds"))

N0_levels <- sort(unique(onesample_power_df$N0))
plt_df <- 
  onesample_power_df %>% 
  filter(method_name == "upstrap") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_N_pointer <- plt_df %>% filter(N0 == N1)
plt_x_axis_max <- plt_df %>% filter(N0 == max(N0)) %>% filter(power_est_aggmedian >= 0.95) %>% pull(N1) %>% min()
plt_x_axis_max

# trick to handle Error: Aesthetics can not vary with a ribbon, so we duplicate
# the bridge observation so as it has both levels
plt_df_dupl <- plt_df %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df)

plt_1 <- 
  ggplot(
    plt_df_dupl, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), 
              alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  geom_line(data = onesample_power_theoret_df, aes(x = N1, y = power_est, group =1), 
            inherit.aes = FALSE, size = 0.4) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  # scale_x_continuous(breaks = c(0, plt_df_N_pointer$N0, 100, 150), limits = c(0, 150)) + 
  scale_x_continuous(breaks = seq(0, 150, by = 50), limits = c(0, 150)) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 1: linear regression with intercept only\n(Est. method: one-sample t-test)") + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp) 
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

plt_df_dupl <- plt_df %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df)

plt_2 <- 
  ggplot(
    plt_df_dupl, 
    aes(x = N1, y = power_est_aggmedian, color = type)
  ) +
  geom_ribbon(aes(ymin = power_est_aggmedian_lwr, ymax = power_est_aggmedian_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  geom_line(data = twosample_power_theoret_df, aes(x = N1, y = power_est, group =1), 
            inherit.aes = FALSE, size = 0.4) + 
  geom_vline(data = plt_df_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_N_pointer, aes(yintercept = power_est_aggmedian), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_N_pointer, aes(x = N1, y = power_est_aggmedian), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  # scale_x_continuous(breaks = c(plt_df_N_pointer$N0, 200, 250, 300)) + 
  scale_x_continuous(breaks = seq(0, 300, by = 50)) + 
  theme(
    legend.position = c(0.8, 0.1),
    legend.title=element_text(size = legend_font_size_numexp),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 2: simple linear regression\n(Est. method: two-sample t-test)",
       color = "Upstrap type:",
       fill = "Upstrap type:"
       ) + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp) 
plt_2


# combine the two plots
plt_list = list(plt_1, plt_2)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

# plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/2021-03-01-ttest_power_est.png")
plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/ttest_power_est.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 12)





