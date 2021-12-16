
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
# our_dir <- paste0(here::here(), "/numerical_experiments/results_figures/2021-08-07")
our_dir <- paste0(here::here(), "/numerical_experiments/results_figures/2021-12-04")
dir.create(our_dir)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 3

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lm_testcoef_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)

# plot output path 
plt_fpath <- paste0(our_dir, "/lm_testcoef_power.png")

eff_tar_levels <- c("0.5", "1",  "observed")
eff_tar_labels <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels <- c("upstrap_power", "comparator_power")
name_labels <- c("upstrap", "SIMR")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 2
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  if (eff_tar_tmp == eff_tar_levels[3]){
    plt_df_tru <- plt_df_tru %>% filter(name == "")
  }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "",
         color = "") + 
    theme(legend.position = "none")  + 
    scale_color_manual(values = my_pal_sub1) 
  ## add legend 
  if (i == length(eff_tar_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 3") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.8)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 4

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glm_testcoef_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)

# plot output path 
plt_fpath <- paste0(our_dir, "/glm_testcoef_power.png")

eff_tar_levels <- c("0.5", "1",  "observed")
eff_tar_labels <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels <- c("upstrap_power", "comparator_power")
name_labels <- c("upstrap", "SIMR")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 2
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  if (eff_tar_tmp == eff_tar_levels[3]){
    plt_df_tru <- plt_df_tru %>% filter(name == "")
  }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "",
         color = "") + 
    theme(legend.position = "none")  + 
    scale_color_manual(values = my_pal_sub1) 
  ## add legend 
  if (i == length(eff_tar_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 4") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.8)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 5

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lmm_testcoef_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)

# plot output path 
plt_fpath <- paste0(our_dir, "/lmm_testcoef_power.png")

eff_tar_levels <- c("0.5", "1",  "observed")
eff_tar_labels <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels <- c("upstrap_power", "comparator_power")
name_labels <- c("upstrap", "SIMR")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 2
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  if (eff_tar_tmp == eff_tar_levels[3]){
    plt_df_tru <- plt_df_tru %>% filter(name == "")
  }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "",
         color = "") + 
    theme(legend.position = "none")  + 
    scale_color_manual(values = my_pal_sub1) 
  ## add legend 
  if (i == length(eff_tar_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 5") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.8)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 6

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glmm_testcoef_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)

# plot output path 
plt_fpath <- paste0(our_dir, "/glmm_testcoef_power.png")

eff_tar_levels <- c("0.5", "1",  "observed")
eff_tar_labels <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels <- c("upstrap_power", "comparator_power")
name_labels <- c("upstrap", "SIMR")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 2
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){
    paste0(N_tar * 2, " (", N_tar, "+", N_tar, ")")
  })
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  if (eff_tar_tmp == eff_tar_levels[3]){
    plt_df_tru <- plt_df_tru %>% filter(name == "")
  }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "",
         color = "") + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + 
    scale_color_manual(values = my_pal_sub1) 
  ## add legend 
  if (i == length(eff_tar_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 6") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.4)

