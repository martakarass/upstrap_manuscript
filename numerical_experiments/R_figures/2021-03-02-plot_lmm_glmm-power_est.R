
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

show_col(my_pal)
show_col(my_pal_numexp)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 3a: GEE linear

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

# trick to handle Error: Aesthetics can not vary with a ribbon, so we duplicate
# the bridge observation so as it has both levels
plt_df_boot_dupl <- plt_df_boot %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df_boot)

plt_3a <- 
  ggplot(
    plt_df_boot_dupl, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # gold standard
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group = 1), size = 0.4, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_lwr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_upr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  # pointer 
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  # scale_x_continuous(breaks = c(0, plt_df_boot_N_pointer$N0, 50, 100, 150)) +
  scale_x_continuous(breaks = seq(0, 250, by = 50), limits = c(NA, 150)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 5: multilevel linear regression\n(Est. method: GEE)",
       color = "Upstrap type:") + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp)
plt_3a

# ------------------------------------------------------------------------------
# PLOT 3b: GEE logistic

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

plt_df_boot_dupl <- plt_df_boot %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df_boot) %>%
  ungroup()
                                
plt_3b <- 
  ggplot(
    plt_df_boot_dupl, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # gold standard
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group = 1), size = 0.4, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_lwr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_upr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  # pointer 
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = seq(0, 250, by = 50), limits = c(NA, 250)) +
  theme(
    legend.position = c(0.8, 0.25),
    legend.title=element_text(size = legend_font_size_numexp),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 6: multilevel generalized linear regression\n(Est. method: GEE)",
       color = "Upstrap type:", fill = "Upstrap type:") + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp)
plt_3b


# combine the two plots
plt_list = list(plt_3a, plt_3b)

plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/lmm_glmm_power_est_GEE.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 5)
rm(plt_3a, plt_3b)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 4a: LMM

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-lmm_trt_agg")
power_df <- readRDS(paste0(fdir_tmp, "/lmm_trt_power_est_bootCI.rds"))

N0_levels <- sort(unique(power_df$N0))
plt_df_boot <- 
  power_df %>% 
  filter(method_name == "out_boot_LMM") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_gold <- 
  power_df %>% 
  filter(method_name == "out_gold_LMM") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_boot_N_pointer <- plt_df_boot %>% filter(N0 == N1)

# trick to handle Error: Aesthetics can not vary with a ribbon, so we duplicate
# the bridge observation so as it has both levels
plt_df_boot_dupl <- plt_df_boot %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df_boot)

plt_4a <- 
  ggplot(
    plt_df_boot_dupl, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # gold standard
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group = 1), size = 0.4, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_lwr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_upr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  # pointer 
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  # scale_x_continuous(breaks = c(0, plt_df_boot_N_pointer$N0, 50, 100, 150)) +
  scale_x_continuous(breaks = seq(0, 250, by = 50), limits = c(NA, 150)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 5: multilevel linear regression\n(Est. method: LMM)",
       color = "Upstrap type:") + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp)
plt_4a

# ------------------------------------------------------------------------------
# PLOT 4b: GLMM

fdir_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-02-17-glmm_trt_agg")
power_df <- readRDS(paste0(fdir_tmp, "/glmm_trt_power_est_bootCI.rds"))

N0_levels <- sort(unique(power_df$N0))
plt_df_boot <- 
  power_df %>% 
  filter(method_name == "out_boot_GLMM") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_gold <- 
  power_df %>% 
  filter(method_name == "out_gold_GLMM") %>%
  mutate(type = ifelse(N1 > N0, "up", "down")) %>%
  mutate(type = factor(type, levels = c("up", "down"))) %>%
  mutate(N0_fct = factor(N0, levels = N0_levels, labels = paste0("Case: observed sample size N=", N0_levels)))
plt_df_boot_N_pointer <- plt_df_boot %>% filter(N0 == N1)

plt_df_boot_dupl <- plt_df_boot %>% group_by(N0) %>% filter(type == "down") %>% filter(N1 == max(N1)) %>% 
  mutate(type = "up") %>% rbind(plt_df_boot) %>%
  ungroup()

plt_4b <- 
  ggplot(
    plt_df_boot_dupl, 
    aes(x = N1, y = est_median, color = type)
  ) +
  geom_ribbon(aes(ymin = est_bootci_median_lwr, ymax = est_bootci_median_upr, fill = type), alpha = 0.2, color = NA) + 
  geom_line(size = 0.7) + 
  # gold standard
  geom_line(data = plt_df_gold, aes(x = N1, y = est_mean, group = 1), size = 0.4, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_lwr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  geom_line(data = plt_df_gold, aes(x = N1, y = est_bootci_mean_upr, group = 1), size = 0.4, linetype = 3, inherit.aes = FALSE) + 
  # pointer 
  geom_vline(data = plt_df_boot_N_pointer, aes(xintercept = N1), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_hline(data = plt_df_boot_N_pointer, aes(yintercept = est_median), alpha = 0.5, size = 0.2, color = "black", linetype = 2) + 
  geom_point(data = plt_df_boot_N_pointer, aes(x = N1, y = est_median), size = 2, color = "black") + 
  facet_wrap(~ N0_fct, ncol = 1) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
  scale_x_continuous(breaks = seq(0, 250, by = 50), limits = c(NA, 250)) +
  theme(
    legend.position = c(0.8, 0.25),
    legend.title=element_text(size = legend_font_size_numexp),
    legend.background = element_rect(linetype="solid", size=0.3),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp)
  ) + 
  labs(x = TeX("Sample size $M_k$"), 
       y = "Estimated power", 
       title = "Problem 6: multilevel generalized linear regression\n(Est. method: GLMM)",
       color = "Upstrap type:", fill = "Upstrap type:") + 
  scale_color_manual(values = my_pal_numexp) + 
  scale_fill_manual(values = my_pal_numexp)
plt_4b


# combine the two plots
plt_list = list(plt_4a, plt_4b)
plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = FALSE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/lmm_glmm_power_est_GLMM.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 5)
