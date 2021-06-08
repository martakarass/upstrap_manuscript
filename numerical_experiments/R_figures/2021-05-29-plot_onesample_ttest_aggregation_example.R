

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)
source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))


# ------------------------------------------------------------------------------
# PLOT 1 

# read simulated example data
results_dir <- paste0(here::here(), "/numerical_experiments/results/2021-05-29-onesample_ttest_aggegating_comparison")
out_df <- readRDS(paste0(results_dir, "/out_df.rds"))

# aggregate data 
out_df_agg <- out_df %>%
  group_by(M_sample_size) %>%
  summarise(
    powerttest_mean = mean(powerttest),
    powerttest_median = median(powerttest),
    upstrap_mean = mean(upstrap),
    upstrap_median = median(upstrap),
    ttest_mean = mean(ttest)
  ) %>%
  pivot_longer(cols = -M_sample_size) %>%
  separate(name, into = c("est_approach", "results_agg_method"))
head(out_df_agg)

# define gold standard power data 
n_arg <- sort(unique(out_df$M_sample_size))
goldstandard_df <- data.frame(
  M_sample_size = n_arg,
  value = power.t.test(n = n_arg, delta = 0.3, sd = 1, type = "one.sample")$power
)

# pivot pre-aggregated data to generate a plot
out_df_L <- 
  out_df %>%
  pivot_longer(cols = -c(rep_idx, M_sample_size), names_to = "est_approach") 
head(out_df_L)

est_approach_levels <- c("powerttest", "upstrap", "ttest")
est_approach_labels <- c("(a) power.t.test", "(b) upstrap", "(c) t.test")

plt_list <- list()
for (i in 1 : length(est_approach_levels)){ # i <- 2
  est_approach_i  <- est_approach_levels[i]
  out_df_agg_i    <- out_df_agg %>% filter(est_approach == est_approach_i)
  out_df_L_i      <- out_df_L %>% filter(est_approach == est_approach_i)
  plt <- ggplot()
  if (i %in% c(1,2)){
    plt <- plt +
      geom_line(data = out_df_L_i %>% filter(rep_idx <= 100),
                aes(x = M_sample_size, y = value, group = rep_idx),
                color = "black", alpha = 0.1, size = 0.3) 
  }
  plt <- plt + 
    geom_line(data = out_df_agg_i, 
              aes(x = M_sample_size, y = value, color = results_agg_method, group = results_agg_method),
              size = 1, inherit.aes = FALSE) + 
    geom_line(data = goldstandard_df, aes(x = M_sample_size, y = value, group = 1), 
              color = "black", inherit.aes = FALSE) + 
    theme(
      legend.title=element_text(size = legend_font_size_numexp),
      plot.title = element_text(size = title_font_size_numexp),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = facetrid_font_size_numexp)
    ) + 
    labs(x = TeX("Sample size $M_k$"),
         y = "",
         color = "",
         title = est_approach_labels[i]) + 
    scale_y_continuous(limits = c(0,1))  + 
    theme(legend.position = c(0.7, 0.3)) 
  # append plot to plots list
  plt_list[[length(plt_list) + 1]] <- plt 
}

plt <- plot_grid(plotlist = plt_list, ncol = 3, align = "v", byrow = TRUE)
plt

plt_fpath <- paste0(here::here(), "/numerical_experiments/results_figures/onesample_ttest_aggegating_comparison_PLOT1.png")
save_plot(filename = plt_fpath, plot = plt, base_width = 10, base_height = 3.5)

