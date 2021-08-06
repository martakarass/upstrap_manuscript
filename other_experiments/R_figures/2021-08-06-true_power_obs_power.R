
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(cowplot)
library(scales)
source(here::here("numerical_experiments/R/config_figures.R"))
source(here::here("numerical_experiments/R/config_utils.R"))

# save plot directory 
our_dir <- paste0(here::here(), "/other_experiments/results_figures/2021-08-06")
dir.create(our_dir)

# read pre-aggregated data 
fpath_tmp <- paste0(here::here(), "/other_experiments/results_CL_shared/2021-08-06-estimate_power.rds")
dat_agg <- readRDS(fpath_tmp)  
nrow(dat_agg)
head(dat_agg)


err_sd_vec <- sort(unique(dat_agg$err_sd))

i <- 3
err_sd_tmp <- err_sd_vec[i]
name_tmp   <- "powerttest"
effsize_tar_tmp <- "observed"

plt_df_est <- dat_agg %>% filter(err_sd == err_sd_tmp, name == name_tmp, effsize_tar == effsize_tar_tmp)
plt_df_tru <- dat_agg %>% filter(err_sd == err_sd_tmp, name == "true", effsize_tar == effsize_tar_tmp)
plt <- 
  ggplot()

