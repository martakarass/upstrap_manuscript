
#' Generate plots showing comparison of upstrap vs. comparator. 

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)
library(ggrepel)
source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 1 

# read simulated raw data -- upstrap
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-onesample_ttest_raw")
fnames_all <- list.files(res_fdir_raw, full.names = TRUE)
fnames_all

N1_min <- 5
N1_max <- 200 
N1_grid <- N1_min : N1_max

exp_N <- 15

# read raw data: upstrap
pat_tmp <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-onesample_ttest_raw/mat_out_upstrap_N0_25.rds")
df_1 <- 
  readRDS(pat_tmp) %>% 
  as.data.frame() %>%
  rename_all(., ~as.character(N1_grid)) %>%
  mutate(method_name = "upstrap") %>%
  mutate(experiment_id = row_number()) %>%
  filter(experiment_id <= exp_N) %>%
  pivot_longer(cols = -c(method_name, experiment_id), names_to = "N1", values_to = "power_est") %>%
  mutate(N1 = as.numeric(N1))
head(df_1)

# read raw data: power.t.test
pat_tmp <- paste0(here::here(), "/numerical_experiments/results_CL/2021-02-17-onesample_ttest_raw/mat_out_powerttest_N0_25.rds")
df_2 <- 
  readRDS(pat_tmp) %>% 
  as.data.frame() %>%
  rename_all(., ~as.character(N1_grid)) %>%
  mutate(method_name = "power.t.test") %>%
  mutate(experiment_id = row_number()) %>%
  filter(experiment_id <= exp_N) %>%
  pivot_longer(cols = -c(method_name, experiment_id), names_to = "N1", values_to = "power_est") %>%
  mutate(N1 = as.numeric(N1))
head(df_2)


plt_list <- list()
N1_val <- 55

plt_tmp <- 
  ggplot(df_1 %>% filter(N1 <= 100), 
         aes(x = N1, y = power_est, group = experiment_id)) + 
  geom_line(aes(color = method_name), alpha = 0.2, size = 0.3) + 
  geom_vline(xintercept = N1_val, size = 0.3, color = "red", linetype = 2) + 
  geom_text_repel(data = df_1 %>% filter(N1 == N1_val),
                  aes(x = N1, y = power_est, label = experiment_id),
                  box.padding = 1, 
                  max.overlaps = Inf, 
                  # segment.curvature = -0.1,
                  # segment.ncp = 2,
                  # segment.angle = 20, 
                  color = "red") +
  theme(
    legend.title=element_text(size = legend_font_size_numexp),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp),
    legend.background = element_rect(size=0.1, linetype="solid", 
                                     colour ="black")
  ) + 
  scale_color_manual(values = "black") + 
  labs(x = TeX("Sample size $M_k$"),
       y = "Estimated power",
       color = "Method: ") + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) + 
  theme(legend.position = c(0.8, 0.3)) 
# plt_tmp
plt_list[[length(plt_list) + 1]] <- plt_tmp 

plt_tmp <- 
  ggplot(df_2 %>% filter(N1 <= 100), 
         aes(x = N1, y = power_est, group = experiment_id)) + 
  geom_line(aes(color = method_name), alpha = 0.2, size = 0.3) + 
  geom_vline(xintercept = N1_val, size = 0.3, color = "red", linetype = 2) + 
  geom_text_repel(data = df_2 %>% filter(N1 == N1_val),
                  aes(x = N1, y = power_est, label = experiment_id),
                  box.padding = 0.8, 
                  max.overlaps = Inf, 
                  # segment.curvature = -0.1,
                  # segment.ncp = 2,
                  # segment.angle = 20, 
                  color = "red") +
  theme(
    legend.title=element_text(size = legend_font_size_numexp),
    plot.title = element_text(size = title_font_size_numexp),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = facetrid_font_size_numexp),
    legend.background = element_rect(size=0.1, linetype="solid", 
                                     colour ="black")
  ) + 
  scale_color_manual(values = "black") + 
  labs(x = TeX("Sample size $M_k$"),
       y = "Estimated power",
       color = "Method: ") + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) + 
  theme(legend.position = c(0.8, 0.3)) 
# plt_tmp
plt_list[[length(plt_list) + 1]] <- plt_tmp 

plt <- plot_grid(plotlist = plt_list, ncol = 2, align = "v", byrow = TRUE)
plot(plt)

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_realizations_upstrap_and_powerttest.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 4)






