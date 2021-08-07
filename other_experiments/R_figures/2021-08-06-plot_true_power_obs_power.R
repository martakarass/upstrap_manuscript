
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)
# source(here::here("numerical_experiments/R/config_figures.R"))
# source(here::here("numerical_experiments/R/config_utils.R"))

# save plot directory 
our_dir <- paste0(here::here(), "/other_experiments/results_figures/2021-08-06")
dir.create(our_dir)

# read pre-aggregated data 
fpath_tmp <- paste0(here::here(), "/other_experiments/results_CL_shared/2021-08-06-estimate_power.rds")
dat_agg <- readRDS(fpath_tmp)  
nrow(dat_agg)
head(dat_agg)

# format
dat_agg <- mutate(dat_agg, N_obs_fct = factor(N_obs))

# params
gg_base_size <- 9
err_sd_vec <- sort(unique(dat_agg$err_sd))

draw_label_text <- c(
  "Compare\n
  (1) 'true power' to detect effect size 0.01,\n
  (2) power to detect observed effect size via power.t.test()",
  "Compare\n
  (1) 'true power' to detect effect size 0.01,\n
  (3) power to detect effect size 0.01 via power.t.test()",
  "Compare\n
  (1) 'true power' to detect effect size 0.01,\n
  (2) power to detect observed effect size via upstrap",
  "Compare\n
  (1) 'true power' to detect effect size 0.01,\n
  (3) power to detect effect size 0.01 via upstrap"
)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 1 

name_tmp   <- "powerttest"
effsize_tar_tmp <- "observed"
N_obs_levels <- c(1000, 2000, 4000)
N_obs_labels <- paste0("n = ", N_obs_levels)

plt_list <- list()
for (i in 1 : 4){
  err_sd_tmp <- err_sd_vec[i]
  # plot data frames
  plt_df_est <- dat_agg %>% filter(err_sd == err_sd_tmp, name == name_tmp, effsize_tar == effsize_tar_tmp)
  plt_df_tru <- dat_agg %>% filter(err_sd == err_sd_tmp, name == "true")
  plt_df_est$N_obs_fct <- factor(plt_df_est$N_obs, levels = N_obs_levels, labels = N_obs_labels) 
  pd <- position_dodge(500)
  plt <- 
    ggplot(plt_df_est, aes(x = N_tar, y = value_q50, color = N_obs_fct)) + 
    geom_point(data = plt_df_tru, aes(x = N_obs, y = value_mean), inherit.aes = FALSE, size = 4) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd) +
    geom_point(size = 2, position = pd) + 
    labs(x = "Target sample size", 
         y = "Power",
         color = "Observed\nsample size",
         title = paste0("SD = ", err_sd_tmp)) + 
    theme_classic(base_size = gg_base_size) + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    scale_x_continuous(breaks = unique(plt_df_est$N_tar))
  # add legend
  if (i == 1){
    plt <- plt +  theme(legend.position = c(0.85, 0.45),
                        legend.title=element_text(size = gg_base_size - 2),
                        legend.background = element_rect(fill = alpha("grey", 0.1)))   
  }
  plt_list[[length(plt_list) + 1]] <- plt
}


plt_all <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt_title <- ggdraw() + 
  draw_label(draw_label_text[1], 
             fontface = 'bold', size = gg_base_size)
# plt_all_F <- plt_all
plt_all_F <- plot_grid(plt_title, plt_all, ncol = 1, rel_heights = c(0.1, 1), scale = 0.95)
plt_all_F

plt_fpath <- paste0(our_dir, "/compare_powerttest_observed.png")
save_plot(filename = plt_fpath, plot = plt_all_F, base_width = 4, base_height  = 9)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 2 

name_tmp   <- "powerttest"
effsize_tar_tmp <- "0.01"
N_obs_levels <- c(1000, 2000, 4000)
N_obs_labels <- paste0("n = ", N_obs_levels)

plt_list <- list()
for (i in 1 : 4){
  err_sd_tmp <- err_sd_vec[i]
  # plot data frames
  plt_df_est <- dat_agg %>% filter(err_sd == err_sd_tmp, name == name_tmp, effsize_tar == effsize_tar_tmp)
  plt_df_tru <- dat_agg %>% filter(err_sd == err_sd_tmp, name == "true")
  plt_df_est$N_obs_fct <- factor(plt_df_est$N_obs, levels = N_obs_levels, labels = N_obs_labels) 
  pd <- position_dodge(500)
  plt <- 
    ggplot(plt_df_est, aes(x = N_tar, y = value_q50, color = N_obs_fct)) + 
    geom_point(data = plt_df_tru, aes(x = N_obs, y = value_mean), inherit.aes = FALSE, size = 4) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd) +
    geom_point(size = 2, position = pd) + 
    labs(x = "Target sample size", 
         y = "Power",
         color = "Observed\nsample size",
         title = paste0("SD = ", err_sd_tmp)) + 
    theme_classic(base_size = gg_base_size) + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    scale_x_continuous(breaks = unique(plt_df_est$N_tar))
  # add legend
  if (i == 1){
    plt <- plt +  theme(legend.position = c(0.85, 0.45),
                        legend.title=element_text(size = gg_base_size - 2),
                        legend.background = element_rect(fill = alpha("grey", 0.1)))   
  }
  plt_list[[length(plt_list) + 1]] <- plt
}


plt_all <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt_title <- ggdraw() + 
  draw_label(draw_label_text[2], 
             fontface = 'bold', size = gg_base_size)
plt_all_F <- plot_grid(plt_title, plt_all, ncol = 1, rel_heights = c(0.1, 1), scale = 0.95) 
# plt_all_F

plt_fpath <- paste0(our_dir, "/compare_powerttest_001.png")
save_plot(filename = plt_fpath, plot = plt_all_F, base_width = 4, base_height  = 9)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 3

name_tmp   <- "upstrap"
effsize_tar_tmp <- "observed"
N_obs_levels <- c(1000, 2000, 4000)
N_obs_labels <- paste0("n = ", N_obs_levels)

plt_list <- list()
for (i in 1 : 4){ # i <- 1
  err_sd_tmp <- err_sd_vec[i]
  # plot data frames
  plt_df_est <- dat_agg %>% filter(err_sd == err_sd_tmp, name == name_tmp, effsize_tar == effsize_tar_tmp)
  plt_df_tru <- dat_agg %>% filter(err_sd == err_sd_tmp, name == "true")
  plt_df_est$N_obs_fct <- factor(plt_df_est$N_obs, levels = N_obs_levels, labels = N_obs_labels) 
  pd <- position_dodge(500)
  plt <- 
    ggplot(plt_df_est, aes(x = N_tar, y = value_q50, color = N_obs_fct)) + 
    geom_point(data = plt_df_tru, aes(x = N_obs, y = value_mean), inherit.aes = FALSE, size = 4) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd) +
    geom_point(size = 2, position = pd) + 
    labs(x = "Target sample size", 
         y = "Power",
         color = "Observed\nsample size",
         title = paste0("SD = ", err_sd_tmp)) + 
    theme_classic(base_size = gg_base_size) + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    scale_x_continuous(breaks = unique(plt_df_est$N_tar))
  # add legend
  if (i == 1){
    plt <- plt +  theme(legend.position = c(0.85, 0.45),
                        legend.title=element_text(size = gg_base_size - 2),
                        legend.background = element_rect(fill = alpha("grey", 0.1)))   
  }
  plt_list[[length(plt_list) + 1]] <- plt
}


plt_all <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt_title <- ggdraw() + 
  draw_label(draw_label_text[3], 
             fontface = 'bold', size = gg_base_size)
plt_all_F <- plot_grid(plt_title, plt_all, ncol = 1, rel_heights = c(0.1, 1), scale = 0.95) 
# plt_all_F

plt_fpath <- paste0(our_dir, "/compare_upstrap_observed.png")
save_plot(filename = plt_fpath, plot = plt_all_F, base_width = 4, base_height  = 9)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 4 

name_tmp   <- "upstrap"
effsize_tar_tmp <- "0.01"
N_obs_levels <- c(1000, 2000, 4000)
N_obs_labels <- paste0("n = ", N_obs_levels)

plt_list <- list()
for (i in 1 : 4){
  err_sd_tmp <- err_sd_vec[i]
  # plot data frames
  plt_df_est <- dat_agg %>% filter(err_sd == err_sd_tmp, name == name_tmp, effsize_tar == effsize_tar_tmp)
  plt_df_tru <- dat_agg %>% filter(err_sd == err_sd_tmp, name == "true")
  plt_df_est$N_obs_fct <- factor(plt_df_est$N_obs, levels = N_obs_levels, labels = N_obs_labels) 
  pd <- position_dodge(500)
  plt <- 
    ggplot(plt_df_est, aes(x = N_tar, y = value_q50, color = N_obs_fct)) + 
    geom_point(data = plt_df_tru, aes(x = N_obs, y = value_mean), inherit.aes = FALSE, size = 4) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd) +
    geom_point(size = 2, position = pd) + 
    labs(x = "Target sample size", 
         y = "Power",
         color = "Observed\nsample size",
         title = paste0("SD = ", err_sd_tmp)) + 
    theme_classic(base_size = gg_base_size) + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    scale_x_continuous(breaks = unique(plt_df_est$N_tar))
  # add legend
  if (i == 1){
    plt <- plt +  theme(legend.position = c(0.85, 0.45),
                        legend.title=element_text(size = gg_base_size - 2),
                        legend.background = element_rect(fill = alpha("grey", 0.1)))   
  }
  plt_list[[length(plt_list) + 1]] <- plt
}


plt_all <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt_title <- ggdraw() + 
  draw_label(draw_label_text[4], 
             fontface = 'bold', size = gg_base_size)
plt_all_F <- plot_grid(plt_title, plt_all, ncol = 1, rel_heights = c(0.1, 1), scale = 0.95) 
# plt_all_F

plt_fpath <- paste0(our_dir, "/compare_upstrap_001.png")
save_plot(filename = plt_fpath, plot = plt_all_F, base_width = 4, base_height  = 9)
