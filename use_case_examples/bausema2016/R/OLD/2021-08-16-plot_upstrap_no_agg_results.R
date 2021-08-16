rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-05-02-upstrap_no_agg_results.rds")
out_df <- readRDS(out_power_df_fpath)
str(out_df)
table(out_df$zone)
table(out_df$time_point)
summary(out_df$B_boot_success)
summary(out_df$eval_time_secs) / (60 * 60) # in hours

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
  out_df_i  <- out_df %>% filter(zone == name_i)
  # if (i == 1) {
  #   out_df_i <- out_df_i %>% filter(power > 0.9)
  # }
  M_min_i <- min(out_df_i %>% filter(power >= 0.8) %>% pull(sample_size_M))
  plt_i <- 
    ggplot(out_df_i %>% filter(sample_size_M < 30), 
           aes(x = sample_size_M, y = power, group = 1)) + 
    geom_hline(yintercept = 0.8, color = "blue", alpha = 0.4) + 
    # annotate('segment', x = M_min_i, y = 0, xend = M_min_i, yend = 0.8, color = "blue", alpha = 0.6, linetype = 2) 
    geom_vline(xintercept = M_min_i, color = "blue", alpha = 0.6, linetype = 2) + 
    geom_point(alpha = 0.8, color = "darkgreen") + 
    geom_line(alpha = 0.8, color = "darkgreen") + 
    # theme(legend.position = "none") + 
    labs(x = TeX("Group size $M_k$"), 
         y = "Power",
         title  = name_labels[i])  + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) + 
    scale_x_continuous(breaks = seq(5, 30, by = 5)) 
  plt_i
  plt_list[[length(plt_list) + 1]] <- plt_i
}
plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = TRUE)
plt

plt_path <- paste0(here::here(), "/use_case_examples/bausema2016/results_figures/2021-05-02-upstrap_no_agg_results.png")
ggsave(filename = plt_path, plot = plt, width = 10, height = 3.3) 