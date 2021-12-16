
#' @description 
#' Generate plot for simulation problems where the focus is to estimate power
#' for situations with a changed proportion of covariate.
#' 
#' Note:
#' - The x-axis ("Target sample size") uses expression "2 * N_obs" which is 
#' equal to total target sample size (both treatment arms); in former simulations,
#' "N_obs" alone was used to denote sample size of each of the two treatment arms. 

rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)
library(latex2exp)

source(here::here("numerical_experiments/R/config_figures_2.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# define color palette based on wider palette defined in config_figures.R 
my_pal_sub1 <- my_pal[c(4,6)]
my_pal_sub2 <- my_pal[c(7,8)]
# show_col(my_pal)

# save plot directory 
our_dir <- paste0(here::here(), "/numerical_experiments/results_figures/2021-12-04")
dir.create(our_dir)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ONE-SAMPLE T-TEST: vary N_obs and error(s) sd 

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-twosample_ttest_vary_covprop_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)
head(dat_agg)

# plot output path 
plt_fpath <- paste0(our_dir, "/twosample_ttest_power_vary_covprop.png")

cov_prop_levels <- c(0.05, 0.10, 0.30)
# cov_prop_labels <- paste0("Target effect\nsize: ", cov_prop_levels)
cov_prop_labels2 <- c(
  TeX("Target proportion of $X_1=1: \\; 0.05$"),
  TeX("Target proportion of $X_1=1: \\; 0.1$"),
  TeX("Target proportion of $X_1=1: \\; 0.3$"))
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")


plt_list <- list()
for (i in 1 : length(cov_prop_levels)){ # i <- 1
  cov_prop_tmp <- cov_prop_levels[i]
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){N_tar * 2}) 
  # get plot data 
  plt_df <- 
    dat_agg %>% 
    filter(cov_prop == cov_prop_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    mutate(N_obs = factor(N_obs)) 
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", cov_prop == cov_prop_tmp) 
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(size = 1.5) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75, color = name), width = .1)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) +
    # scale_x_continuous(breaks = x_axis_breaks) + 
    annotate("text", x = 5, y = 1.1, label = cov_prop_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") + 
    theme(legend.position = "none") +
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + 
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 7") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  if (i == length(cov_prop_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.8)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LM: vary N_obs and error(s) sd 

# data input directory 
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lm_testcoef_vary_covprop_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
table(dat_agg$name)
head(dat_agg)

# plot output path 
plt_fpath <- paste0(our_dir, "/lm_testcoef_power_vary_covprop.png")


cov_prop_levels <- c(0.05, 0.10, 0.30)
# cov_prop_labels <- paste0("Target effect\nsize: ", cov_prop_levels)
cov_prop_labels2 <- c(
  TeX("Target proportion of $X_1=1: \\; 0.05$"),
  TeX("Target proportion of $X_1=1: \\; 0.1$"),
  TeX("Target proportion of $X_1=1: \\; 0.3$"))
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")


plt_list <- list()
for (i in 1 : length(cov_prop_levels)){ # i <- 1
  cov_prop_tmp <- cov_prop_levels[i]
  # get x-axis breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){
    N_tar * 2
  })  # get plot data 
  plt_df <- 
    dat_agg %>% 
    filter(cov_prop == cov_prop_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    mutate(N_obs = factor(N_obs)) 
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", cov_prop == cov_prop_tmp) 
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(size = 1.5) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75, color = name), width = .1)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) +
    # scale_x_continuous(breaks = x_axis_breaks) + 
    annotate("text", x = 5, y = 1.1, label = cov_prop_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") + 
    theme(legend.position = "none") +
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + 
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 8") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  if (i == length(cov_prop_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 2.8)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LMM: vary N_obs and error(s) sd 

# data input directory
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lmm_testcoef_vary_covprop_agg.rds")
dat_agg <- readRDS(fpath_tmp)
table(dat_agg$name)
head(dat_agg)

# plot output path
plt_fpath <- paste0(our_dir, "/lmm_testcoef_power_vary_covprop.png")

sd_sigma_levels <- sort(unique(dat_agg$sd_sigma))
sd_sigma_labels2 <- paste0("Error term std: ", sd_sigma_levels)
name_levels <- c("upstrap_power")
name_labels <- c("upstrap")


plt_list <- list()
for (i in 1 : length(cov_prop_levels)){ # i <- 1
  cov_prop_tmp <- cov_prop_levels[i]
  # get x-axos breaks
  x_axis_breaks <- unique(dat_agg$N_tar)
  x_axis_labels <- sapply(x_axis_breaks, function(N_tar){
    size_0 <- round(N_tar * 2 * (1-cov_prop_tmp))
    size_1 <- (N_tar * 2) - size_0
    paste0(N_tar * 2, " (", size_0, "+", size_1, ")")
  })
  pd <- position_dodge(25)
  # get plot data 
  plt_df <- 
    dat_agg %>% 
    filter(cov_prop == cov_prop_tmp, name %in% name_levels) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels)) %>%
    mutate(N_obs = factor(N_obs)) 
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", cov_prop == cov_prop_tmp) 
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(size = 1.5) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75, color = name), width = .1)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels) +
    # scale_x_continuous(breaks = x_axis_breaks) + 
    annotate("text", x = 5, y = 1.1, label = cov_prop_labels2[i], 
             hjust = 0, size = 4) + 
    labs(x = TeX("Target sample size"),
         y = "Power",
         shape = "Observed\nsample size") + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  +
    scale_color_manual(values = my_pal_sub1)
  if (i == 1){
    plt <- plt + labs(title = "Simulation problem 9") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  if (i == length(cov_prop_levels)){
    plt <- plt + 
      theme(legend.position = c(0.78, 0.2),
            legend.title=element_blank(),
            legend.background = element_rect(fill=alpha('white', 0.6), color = NA))  
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = FALSE)
# plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.4)

