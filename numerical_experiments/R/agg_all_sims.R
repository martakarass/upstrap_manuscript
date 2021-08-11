
#' This script aggregate power of rejecting H0 in all the problems. 

rm(list = ls())
library(here)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


#' Function to read and aggregate data for simulation problems 1,2
read_and_agg <- function(fdir_raw){
  # read and combine data 
  fnames <- list.files(fdir_raw, full.names = TRUE)
  fnames <- fnames[grepl("arrayjob", fnames)]
  dat_l  <- lapply(fnames, readRDS)
  dat    <- do.call("rbind", dat_l)
  message(paste0("Read raw files from dir: ", fdir_raw))
  message(paste0("Raw files count: ", length(fnames)))
  # mutate the data 
  dat    <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))
  dat    <- mutate(dat, name = recode(name, powerttest_power = "comparator_power"))
  dat    <- mutate(dat, name = recode(name, simr_power = "comparator_power"))
  # aggregate data: power mean 
  dat_agg <- 
    dat %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = sum(!is.na(value)),
      value_mean = mean(value),
      value_sd = sd(value),
      value_q25 = quantile(value, 0.25),
      value_q50 = median(value),
      value_q75 = quantile(value, 0.75)
    ) %>%
    ungroup()
  # aggregate data: power errors 
  dat_truepower_agg <- 
    dat %>%
    filter(name == "true_power") %>% 
    group_by(N_tar, N_obs, name, eff_tar) %>%
    summarise(value = mean(value)) %>% 
    pivot_wider(names_from = name, values_from = value) %>%
    ungroup() %>%
    select(-N_obs)
  dat_diff <- 
    dat %>%
    filter(name != "true_power", eff_tar != "observed") %>% 
    pivot_wider(names_from = name, values_from = value) %>%
    left_join(dat_truepower_agg, by = c("N_tar", "eff_tar")) %>%
    mutate(
      # upstrap
      pe_ups_true =  100 * (upstrap_power - true_power)/true_power,
      pe_cmp_true =  100 * (comparator_power - true_power)/true_power,
      pe_ups_cmp  =  100 * (upstrap_power - comparator_power)/comparator_power
      ) %>%
    select(-c(upstrap_power, comparator_power, true_power)) %>%
    pivot_longer(cols = starts_with("pe"))
  dat_agg_diff <- 
    dat_diff %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = sum(!is.na(value)),
      value_mean = mean(value),
      value_sd = sd(value),
      value_q25 = quantile(value, 0.25),
      value_q50 = median(value),
      value_q75 = quantile(value, 0.75)
    ) %>%
    ungroup()
  # rbind 
  dat_out <- dat_agg %>% rbind(dat_agg_diff) 
  # print messages
  message(paste0("Agg file nrows: ", nrow(dat_out)))
  dat_out
}


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 1 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-onesample_ttest_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-onesample_ttest_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
# head(out_tmp)
# table(out_tmp$name)

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 2 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-twosample_ttest_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-twosample_ttest_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
# head(out_tmp)
# out_tmp %>% filter(is.na(value_mean))

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 3 -- lm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-lm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lm_testcoef_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
# head(out_tmp)
# table(out_tmp$name)

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 4 -- glm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-glm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glm_testcoef_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
# head(out_tmp)

# Raw files count: 1000

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 5 -- lmm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-lmm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-lmm_testcoef_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
# head(out_tmp)

# Raw files count: 811

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 6 -- glmm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-07-glmm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-07-glmm_testcoef_agg.rds")
out_tmp   <- read_and_agg(fdir_raw)
head(out_tmp) 

# Raw files count: 390

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)

#' cd $ups
#' git add --all
#' git commit -m 'add updated results from the cluster'
#' git push


