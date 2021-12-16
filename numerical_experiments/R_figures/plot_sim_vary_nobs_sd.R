
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures_2.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# define color palette based on wider palette defined in config_figures.R 
my_pal_sub1 <- my_pal[c(4,6)]
my_pal_sub2 <- my_pal[c(7,8)]
# show_col(my_pal)
# show_col(my_pal_sub1)
# show_col(my_pal_sub2)

# save plot directory 
our_dir <- paste0(here::here(), "/numerical_experiments/results_figures/2021-12-04")
dir.create(our_dir)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ONE-SAMPLE T-TEST: vary N_obs and error(s) sd 

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-twosample_ttest_vary_nobs_sd_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)
head(dat_agg)

# plot output path 
plt_fpath <- paste0(our_dir, "/twosample_ttest_power_vary_nobs_sd.png")

sd_sigma_levels <- sort(unique(dat_agg$sd_sigma))
# sd_sigma_labels <- paste0("Target effect\nsize: ", sd_sigma_levels)
sd_sigma_labels2 <- paste0("Error term std: ", sd_sigma_levels)
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")


plt_list <- list()
for (i in 1 : length(sd_sigma_levels)){ # i <- 1
  sd_sigma_tmp <- sd_sigma_levels[i]
  pd <- position_dodge(25)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  N_obs_levels = sort(unique(dat_agg$N_obs))
  N_obs_labels <- sapply(N_obs_levels, function(N_obs){
    paste0(N_obs * 2, " (", N_obs, "+", N_obs, ")")
  })  
  plt_df <- 
    dat_agg %>% 
    filter(sd_sigma == sd_sigma_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    mutate(N_obs = factor(N_obs, levels = N_obs_levels, labels = N_obs_labels)) 
    # mutate(N_obs = factor(N_obs)) 
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", sd_sigma == sd_sigma_tmp) 
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, shape = N_obs)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 1.5, color = my_pal_sub1[1]) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd, color = my_pal_sub1[1])  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    annotate("text", x = 5, y = 1.15, label = sd_sigma_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") + 
    theme(legend.position = "none")  + 
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt +
      theme(legend.position = c(0.75, 0.42),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.background = element_rect(fill=alpha('white', 0.8), color = NA))
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 10") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.9)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LM: vary N_obs and error(s) sd 

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lm_testcoef_vary_nobs_sd_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)
head(dat_agg)

# plot output path 
plt_fpath <- paste0(our_dir, "/lm_testcoef_power_vary_nobs_sd.png")

sd_sigma_levels <- sort(unique(dat_agg$sd_sigma))
sd_sigma_labels2 <- paste0("Error term std: ", sd_sigma_levels)
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")

plt_list <- list()
for (i in 1 : length(sd_sigma_levels)){ # i <- 1
  sd_sigma_tmp <- sd_sigma_levels[i]
  pd <- position_dodge(25)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  N_obs_levels = sort(unique(dat_agg$N_obs))
  N_obs_labels <- sapply(N_obs_levels, function(N_obs){
    paste0(N_obs * 2, " (", N_obs, "+", N_obs, ")")
  })  
  plt_df <- 
    dat_agg %>% 
    filter(sd_sigma == sd_sigma_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    mutate(N_obs = factor(N_obs)) 
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", sd_sigma == sd_sigma_tmp) 
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, shape = N_obs)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 1.5, color = my_pal_sub1[1]) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd, color = my_pal_sub1[1])  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    annotate("text", x = 5, y = 1.15, label = sd_sigma_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") + 
    theme(legend.position = "none")  + 
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt +
      theme(legend.position = c(0.75, 0.42),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.background = element_rect(fill=alpha('white', 0.8), color = NA))
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 11") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.9)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LMM: vary N_obs and error(s) sd 

# data input directory
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lmm_testcoef_vary_nobs_sd_agg.rds")
dat_agg <- readRDS(fpath_tmp)
table(dat_agg$name)
head(dat_agg)

# plot output path
plt_fpath <- paste0(our_dir, "/lmm_testcoef_power_vary_nobs_sd.png")

sd_sigma_levels <- sort(unique(dat_agg$sd_sigma))
sd_sigma_labels2 <- paste0("Error term std: ", sd_sigma_levels)
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")

plt_list <- list()
for (i in 1 : length(sd_sigma_levels)){ # i <- 1
  sd_sigma_tmp <- sd_sigma_levels[i]
  pd <- position_dodge(25)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){
    paste0(N_tar * 2, " (", N_tar, "+", N_tar, ")")
  })
  N_obs_levels = sort(unique(dat_agg$N_obs))
  N_obs_labels <- sapply(N_obs_levels, function(N_obs){
    paste0(N_obs * 2, " (", N_obs, "+", N_obs, ")")
  })  
  plt_df <-
    dat_agg %>%
    filter(sd_sigma == sd_sigma_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    # mutate(N_obs = factor(N_obs))  %>%
    mutate(N_obs = factor(N_obs, levels = N_obs_levels, labels = N_obs_labels)) 
  plt_df_tru <-
    dat_agg %>%
    filter(name == "true_power", sd_sigma == sd_sigma_tmp)
  plt <-
    ggplot(plt_df, aes(x = N_tar, y = value_q50, shape = N_obs)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE,
               size = 3, shape = 1) +
    geom_point(position = pd, size = 1.5, color = my_pal_sub1[1]) +
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd, color = my_pal_sub1[1])  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) +
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    annotate("text", x = 5, y = 1.15, label = sd_sigma_labels2[i],
             hjust = 0, size = 4) +
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + 
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt +
      theme(legend.position = c(0.75, 0.42),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.background = element_rect(fill=alpha('white', 0.8), color = NA))
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 12") +
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.45)

