
#' Script to aggregate the data

rm(list = ls())
library(tidyverse)

# out_dir <- paste0(here::here(), "/other_experiments/results_CL_shared/2021-08-05-estimate_power")
out_dir <- paste0(here::here(), "/other_experiments/results_CL/2021-08-06-estimate_power")
fnames <- list.files(out_dir, full.names = TRUE)
length(fnames)

dat_l <- lapply(fnames, function(fname){
  dat <- readRDS(fname)
  dat$basename <- basename(fname)
  return(dat)
}) 
dat <- do.call("rbind", dat_l)
length(unique(dat$basename))

dat_agg <- 
  dat %>% 
  mutate(effsize_tar = ifelse(is.na(effsize_tar), "observed", effsize_tar)) %>%
  group_by(name, err_sd, N_tar, N_obs, effsize_tru, effsize_tar) %>%
  summarise(
    cnt = n(),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75)
  ) %>%
  arrange(name, err_sd, N_tar, N_obs, effsize_tru, effsize_tar) %>%
  as.data.frame()
head(dat_agg)

# save results to file 
saveRDS(object = dat_agg, file =  paste0(here::here(), "/other_experiments/results_CL_shared/2021-08-06-estimate_power.rds"))



