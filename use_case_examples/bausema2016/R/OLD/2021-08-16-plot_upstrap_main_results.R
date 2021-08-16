
#' For formatting the data, this script reuses code attached to the following paper: 
#' 
#' - Paper citation: Naudet F, Sakarovitch C, Janiaud P, Cristea I, Fanelli D, Moher D et al. Data sharing and reanalysis of randomized controlled trials in leading biomedical journals with a full data sharing policy: survey of studies published in The BMJ and PLOS Medicine BMJ 2018; 360 :k400 doi:10.1136/bmj.k400
#' - Paper URL: https://www.bmj.com/content/360/bmj.k400
#' - Code URL: https://osf.io/7ghfa/

rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-04-22-upstrap_main_results.rds")
out_df <- readRDS(out_power_df_fpath)
str(out_df)
table(out_df$zone)

out_df %>% filter(sample_size_M == 5) %>% pull(power) %>% round(2)

name_levels <- c("hotspot", "eval_zone_1", "eval_zone_2")
name_labels <- c("Hotspot", "Eval. zone 1-249 m", "Eval. zone 250-500 m")

theme_ggpr <- function(){ 
  font <- "Arial"  
  theme_minimal(base_size = 12) %+replace%    
    theme(panel.grid.major = element_line(size = 0.3),  
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=12))}
theme_set(theme_ggpr())

plt_list <- list()
for (i in 1 : length(name_levels)){ # i <- 1
  name_i  <- name_levels[i]
  out_df_i    <- out_df %>% filter(zone == name_i)
  M_min_i <- min(out_df_i %>% filter(power >= 0.8) %>% pull(sample_size_M))
  plt_i <- 
    ggplot(out_df_i %>% filter(sample_size_M <= 30), 
           aes(x = sample_size_M, y = power, group = 1)) + 
    geom_hline(yintercept = 0.8, color = "blue", alpha = 0.4) + 
    # annotate('segment', x = M_min_i, y = 0, xend = M_min_i, yend = 0.8, color = "blue", alpha = 0.6, linetype = 2) 
    geom_vline(xintercept = M_min_i, color = "blue", alpha = 0.6, linetype = 2) + 
    geom_point(alpha = 0.6) + 
    geom_line(alpha = 0.6) + 
    # theme(legend.position = "none") + 
    labs(x = "Target sample size", 
         y = "Power",
         title  = name_labels[i])  + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) + 
    scale_x_continuous(breaks = seq(5, 30, by = 5)) 
  # plt_i
  plt_list[[length(plt_list) + 1]] <- plt_i
}
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = TRUE)
plt

# plt_path <- paste0(here::here(), "/use_case_examples/bausema2016/results_figures/fig1.png")
plt_path <- paste0(here::here(), "/use_case_examples/bausema2016/results_figures/upstrap_main_results.png")
ggsave(filename = plt_path, plot = plt, width = 10, height = 3.3) 

