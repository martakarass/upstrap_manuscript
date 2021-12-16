
#' For formatting the data, this script reuses code attached to the following paper: 
#' 
#' - Paper citation: Naudet F, Sakarovitch C, Janiaud P, Cristea I, Fanelli D, Moher D et al. Data sharing and reanalysis of randomized controlled trials in leading biomedical journals with a full data sharing policy: survey of studies published in The BMJ and PLOS Medicine BMJ 2018; 360 :k400 doi:10.1136/bmj.k400
#' - Paper URL: https://www.bmj.com/content/360/bmj.k400
#' - Code URL: https://osf.io/7ghfa/

rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 1 

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-04-22-upstrap_main_results.rds")
out_df <- readRDS(out_power_df_fpath) 
# fix 2021-12-11 to have "sample size" denote all units (not: units per treatment arm)
out_df <- mutate(out_df, sample_size_M = 2 * sample_size_M) 
str(out_df)
table(out_df$zone)

out_df %>% filter(sample_size_M == 10) %>% mutate(power = round(power, 2)) 
# out_df %>% filter(sample_size_M == 5) %>% pull(power) %>% round(2)

name_levels <- c("hotspot", "eval_zone_1", "eval_zone_2")
name_labels <- c("Hotspot", "Eval. zone 1-249 m", "Eval. zone 250-500 m")

theme_ggpr <- function(){ 
  font <- "Arial"  
    theme_bw(base_size = 14) %+replace%    
    theme(legend.background = element_rect(fill=alpha('white', 0.6), color = NA),
          panel.grid.major = element_line(size = 0.3),  
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) 
}
theme_set(theme_ggpr())

plt_list <- list()
for (i in 1 : length(name_levels)){ # i <- 1
  name_i  <- name_levels[i]
  out_df_i    <- out_df %>% filter(zone == name_i)
  M_min_i <- min(out_df_i %>% filter(power >= 0.8) %>% pull(sample_size_M))
  plt_i <- 
    ggplot(out_df_i %>% filter(sample_size_M <= 60), 
           aes(x = sample_size_M, y = power, group = 1)) + 
    geom_hline(yintercept = 0.8, color = "blue", alpha = 0.4) + 
    geom_segment(x = M_min_i, xend = M_min_i, y = 0, yend = 1, color = "blue", alpha = 0.6, linetype = 2, inherit.aes = FALSE) +
    geom_point(alpha = 0.6) + 
    geom_line(alpha = 0.6) + 
    # theme(legend.position = "none") + 
    labs(x = "Target sample size", 
         y = "Power",
         title  = "")  + 
    annotate("text", x = 5, y = 1.1, label = name_labels[i], 
             hjust = 0, size = 4) + 
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, by = 0.2)) + 
    scale_x_continuous(breaks = seq(10, 60, by = 10)) 
  if (i == 1){
    plt_i <- 
      plt_i + labs(title = "T-test power") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt_i
}
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = TRUE)
plt

plt_path <- paste0(here::here(), "/use_case_examples/bausema2016/results_figures/upstrap_rda_bousema_ttest.png")
ggsave(filename = plt_path, plot = plt, width = 10, height = 3.3) 



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT 2

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-05-02-upstrap_no_agg_results.rds")
out_df <- readRDS(out_power_df_fpath)
# fix 2021-12-11 to have "sample size" denote all units (not: units per treatment arm)
out_df <- mutate(out_df, sample_size_M = 2 * sample_size_M) 

plt_list <- list()
for (i in 1 : length(name_levels)){ # i <- 1
  name_i  <- name_levels[i]
  out_df_i  <- out_df %>% filter(zone == name_i)
  M_min_i <- min(out_df_i %>% filter(power >= 0.8) %>% pull(sample_size_M))
  plt_i <- 
    ggplot(out_df_i %>% filter(sample_size_M <= 60), 
           aes(x = sample_size_M, y = power, group = 1)) + 
    geom_hline(yintercept = 0.8, color = "blue", alpha = 0.4) + 
    geom_segment(x = M_min_i, xend = M_min_i, y = 0, yend = 1, color = "blue", alpha = 0.6, linetype = 2, inherit.aes = FALSE) +
    geom_point(alpha = 0.6) + 
    geom_line(alpha = 0.6) + 
    # theme(legend.position = "none") + 
    labs(x = "Target sample size", 
         y = "Power",
         title  = "")  +
    annotate("text", x = 5, y = 1.1, label = name_labels[i], 
             hjust = 0, size = 4) + 
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, by = 0.2)) + 
    scale_x_continuous(breaks = seq(10, 60, by = 10)) 
  if (i == 1){
    plt_i <- 
      plt_i + labs(title = "GEE test power") + 
      theme(plot.title = element_text(size = 14, face = "bold"))
  }
  plt_list[[length(plt_list) + 1]] <- plt_i
}
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "hv", byrow = TRUE)
plt

plt_path <- paste0(here::here(), "/use_case_examples/bausema2016/results_figures/upstrap_rda_bousema_gee.png")
ggsave(filename = plt_path, plot = plt, width = 10, height = 3.3) 


