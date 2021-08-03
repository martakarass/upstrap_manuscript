
#' This script aggregate power of rejecting H0 in all the problems. 

rm(list = ls())
library(here)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Function to read and aggregate data for simulation problems 1,2
read_and_agg_A <- function(fdir_raw){
  # read and combine data 
  fnames <- list.files(fdir_raw, full.names = TRUE)
  fnames <- fnames[grepl("arrayjob", fnames)]
  dat_l  <- lapply(fnames, readRDS)
  dat    <- do.call("rbind", dat_l)
  message(paste0("Read raw files from dir: ", fdir_raw))
  message(paste0("Raw files count: ", length(fnames)))
  # mutate the data 
  dat    <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))
  # aggregate data: power mean 
  dat_agg <- 
    dat %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = n(),
      value_mean = mean(value),
      value_sd = sd(value)
    ) %>%
    ungroup()
  # aggregate data: power errors 
  dat_diff <- 
    dat %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(
      ups_cmp_power_diff = upstrap_power - powerttest_power,
      ups_cmp_power_pe = 100 * (upstrap_power - powerttest_power)/powerttest_power,
      ups_cmp_power_ape = abs(ups_cmp_power_pe)) %>%
    select(-c(upstrap_power, powerttest_power)) %>%
    pivot_longer(cols = starts_with("ups_cmp_power"))
  dat_agg_diff <- 
    dat_diff %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = n(),
      value_mean = mean(value),
      value_sd = sd(value)
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

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-onesample_ttest_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-onesample_ttest_agg.rds")
out_tmp   <- read_and_agg_A(fdir_raw)
head(out_tmp)

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 2 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-twosample_ttest_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-twosample_ttest_agg.rds")
out_tmp   <- read_and_agg_A(fdir_raw)
head(out_tmp)

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Function to read and aggregate data for simulation problems 3-6
read_and_agg_B <- function(fdir_raw){
  # read and combine data 
  fnames <- list.files(fdir_raw, full.names = TRUE)
  fnames <- fnames[grepl("arrayjob", fnames)]
  dat_l  <- lapply(fnames, readRDS)
  dat    <- do.call("rbind", dat_l)
  message(paste0("Read raw files from dir: ", fdir_raw))
  message(paste0("Raw files count: ", length(fnames)))
  if (length(fnames) == 0) return(NULL)
  # mutate the data 
  dat    <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))
  dat    <- filter(dat, name != "indepsample_power")
  # aggregate data: power mean 
  dat_agg <- 
    dat %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = n(),
      value_mean = mean(value),
      value_sd = sd(value)
    ) %>%
    ungroup()
  # aggregate data: power errors
  dat_diff <- 
    dat %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(
      ups_cmp_power_diff = upstrap_power - simr_power,
      ups_cmp_power_pe = 100 * (upstrap_power - simr_power)/simr_power,
      ups_cmp_power_ape = abs(ups_cmp_power_pe)) %>%
    select(-c(upstrap_power, simr_power)) %>%
    pivot_longer(cols = starts_with("ups_cmp_power"))
  dat_agg_diff <- 
    dat_diff %>%
    group_by(N_tar, N_obs, name, eff_tru, eff_tar) %>%
    summarise(
      cnt = n(),
      value_mean = mean(value),
      value_sd = sd(value)
    ) %>%
    ungroup()
  dat_out <- rbind(dat_agg, dat_agg_diff)
  message(paste0("Agg file nrows: ", nrow(dat_out)))
  return(dat_out)
}


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 3 -- lm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-07-lm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-07-lm_testcoef_agg.rds")
out_tmp   <- read_and_agg_B(fdir_raw)
head(out_tmp)

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 4 -- glm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-02-glm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-02-glm_testcoef_agg.rds")
out_tmp   <- read_and_agg_B(fdir_raw)
head(out_tmp)
Sys.time()
# "2021-08-02 18:53:11 EDT"; Raw files count: 122

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 5 -- lmm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-07-08-lmm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-07-08-lmm_testcoef_agg.rds")
out_tmp   <- read_and_agg_B(fdir_raw)
head(out_tmp)
Sys.time()
# "2021-08-02 18:53:11 EDT"; Raw files count: 450

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)


# ------------------------------------------------------------------------------
# SIMULATION PROBLEM 6 -- glmm

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-08-02-glmm_testcoef_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-08-02-glmm_testcoef_agg.rds")
out_tmp   <- read_and_agg_B(fdir_raw)
head(out_tmp)
Sys.time()

# save to file
saveRDS(out_tmp, fpath_out)

# remove files 
rm(fdir_raw, fpath_out, out_tmp)

# git add --all
# git commit -m 'add updated results from the cluster'
# git push






