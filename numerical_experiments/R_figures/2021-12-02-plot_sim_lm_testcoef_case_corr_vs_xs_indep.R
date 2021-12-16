
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)

source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# define color palette based on wider palette defined in config_figures.R 
my_pal_sub1 <- my_pal[c(4,6)]
my_pal_sub2 <- my_pal[c(7,8)]
# show_col(my_pal)
# show_col(my_pal_sub1)
# show_col(my_pal_sub2)

# save plot directory 
our_dir <- paste0(here::here(), "/numerical_experiments/results_figures/2021-12-02")
dir.create(our_dir)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION -- CASE (1): independent Xs

# data input  
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-02-lm_testcoef_case_xs_indep_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
head(dat_agg)

# data output
plt_fpath <- paste0(our_dir, "/lm_testcoef_case_xs_indep_power.png")

eff_tar_levels  <- c("0.5")
eff_tar_labels  <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels     <- c("upstrap_power")
name_labels     <- c("upstrap")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 1
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  # if (eff_tar_tmp == eff_tar_levels[3]){
  #   plt_df_tru <- plt_df_tru %>% filter(name == "")
  # }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = unique(plt_df$N_tar)) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4, fill = "green") + 
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
    plt <- plt + labs(title = "Simulation problem: independent Xs") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.0)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SIMULATION -- CASE (2): correlated Xs

# data input  
fpath_tmp <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-02-lm_testcoef_case_corr_agg.rds")
dat_agg <- readRDS(fpath_tmp)  
head(dat_agg)

# data output
plt_fpath <- paste0(our_dir, "/lm_testcoef_case_corr_power.png")

eff_tar_levels  <- c("0.5")
eff_tar_labels  <- paste0("Target effect\nsize: ", eff_tar_levels)
eff_tar_labels2 <- paste0("Target effect size: ", eff_tar_levels)
name_levels     <- c("upstrap_power")
name_labels     <- c("upstrap")
N_obs_tmp <- 50

plt_list <- list()
for (i in 1 : length(eff_tar_levels)){ # i <- 1
  eff_tar_tmp <- eff_tar_levels[i]
  pd <- position_dodge(10)
  plt_df <- 
    dat_agg %>% 
    filter(eff_tar == eff_tar_tmp, name %in% name_levels, N_obs == N_obs_tmp) %>%
    mutate(name = factor(name, levels = name_levels, labels = name_labels))
  plt_df_tru <- 
    dat_agg %>% 
    filter(name == "true_power", eff_tar == eff_tar_tmp) 
  # if (eff_tar_tmp == eff_tar_levels[3]){
  #   plt_df_tru <- plt_df_tru %>% filter(name == "")
  # }
  plt <- 
    ggplot(plt_df, aes(x = N_tar, y = value_q50, color = name)) +
    geom_point(data = plt_df_tru, aes(x = N_tar, y = value_mean), inherit.aes = FALSE, 
               size = 3, shape = 1) + 
    geom_point(position = pd, size = 2) + 
    geom_errorbar(aes(ymin = value_q25, ymax = value_q75), width = .1, position = pd)  +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, NA)) + 
    scale_x_continuous(breaks = unique(plt_df$N_tar)) + 
    # coord_cartesian(ylim = c(0,1)) + 
    annotate("text", x = 5, y = 1.1, label = eff_tar_labels2[i], 
             hjust = 0, size = 4, fill = "green") + 
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
    plt <- plt + labs(title = "Simulation problem: correlated Xs") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt
}

plt <- plot_grid(plotlist = plt_list, ncol = 1, align = "hv", byrow = FALSE)
plt

save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.0)

