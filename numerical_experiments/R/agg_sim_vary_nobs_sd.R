
#' This script aggregate power of rejecting H0 in all the problems. 

rm(list = ls())
library(here)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TWO-SAMPLE T-TEST: vary N_obs and error(s) sd 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-twosample_ttest_vary_nobs_sd_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-twosample_ttest_vary_nobs_sd_agg.rds")

# read and combine data 
fnames <- list.files(fdir_raw, full.names = TRUE)
fnames <- fnames[grepl("arrayjob", fnames)]
message(paste0("Raw files count: ", length(fnames)))
dat_l  <- lapply(fnames, readRDS)
dat    <- do.call("rbind", dat_l)
dim(dat)
head(dat)
# mutate the data 
# dat    <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))
# dat    <- mutate(dat, name = recode(name, powerttest_power = "comparator_power"))
# dat    <- mutate(dat, name = recode(name, simr_power = "comparator_power"))
# aggregate data: estimated power  
dat_agg <- 
  dat %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75),
    prop_success_mean = mean(prop_success),
    prop_success_min = min(prop_success),
    prop_success_max = max(prop_success)
  ) %>%
  ungroup()
summary(dat_agg$prop_success_mean)
summary(dat_agg$prop_success_min)
dat_agg <- dat_agg %>% select(-c(prop_success_mean, prop_success_min, prop_success_max))
dim(dat_agg)
# aggregate data: true mean 
dat_truepower_agg <- 
  dat %>%
  filter(name == "true_power") %>% 
  group_by(N_tar, N_obs, name, eff_tar, sd_sigma) %>%
  summarise(value = mean(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  ungroup() %>%
  select(-N_obs)
# compare with true power
dat_diff_A <- 
  dat %>%
  filter(name != "true_power", eff_tar != "observed") %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(dat_truepower_agg, by = c("N_tar", "eff_tar", "sd_sigma")) %>%
  mutate(
    pe_ups_true =  100 * (upstrap_power - true_power)/true_power
  ) %>%
  select(-c(upstrap_power, true_power)) %>%
  pivot_longer(cols = starts_with("pe"))
dim(dat_diff_A)
dat_diff <- dat_diff_A
dat_agg_diff <- 
  dat_diff %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75)
  ) %>%
  ungroup()
# combine
dat_out <- dat_agg %>% rbind(dat_agg_diff) 

# save 
saveRDS(dat_out, fpath_out)
rm(dat_out)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LM: vary N_obs and error(s) sd 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-lm_testcoef_vary_nobs_sd_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lm_testcoef_vary_nobs_sd_agg.rds")

# read and combine data 
fnames <- list.files(fdir_raw, full.names = TRUE)
fnames <- fnames[grepl("arrayjob", fnames)]
message(paste0("Raw files count: ", length(fnames)))
dat_l  <- lapply(fnames, readRDS)
dat    <- do.call("rbind", dat_l)
dim(dat)
head(dat)
# mutate the data 
# dat    <- mutate(dat, eff_tar = ifelse(is.na(eff_tar), "observed", eff_tar))
# dat    <- mutate(dat, name = recode(name, powerttest_power = "comparator_power"))
# dat    <- mutate(dat, name = recode(name, simr_power = "comparator_power"))
# aggregate data: estimated power  
dat_agg <- 
  dat %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75),
    prop_success_mean = mean(prop_success),
    prop_success_min = min(prop_success),
    prop_success_max = max(prop_success)
  ) %>%
  ungroup()
summary(dat_agg$prop_success_mean)
summary(dat_agg$prop_success_min)
dat_agg <- dat_agg %>% select(-c(prop_success_mean, prop_success_min, prop_success_max))
dim(dat_agg)
# aggregate data: true mean 
dat_truepower_agg <- 
  dat %>%
  filter(name == "true_power") %>% 
  group_by(N_tar, N_obs, name, eff_tar, sd_sigma) %>%
  summarise(value = mean(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  ungroup() %>%
  select(-N_obs)
# compare with true power
dat_diff_A <- 
  dat %>%
  filter(name != "true_power", eff_tar != "observed") %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(dat_truepower_agg, by = c("N_tar", "eff_tar", "sd_sigma")) %>%
  mutate(
    pe_ups_true =  100 * (upstrap_power - true_power)/true_power
  ) %>%
  select(-c(upstrap_power, true_power)) %>%
  pivot_longer(cols = starts_with("pe"))
dim(dat_diff_A)
dat_diff <- dat_diff_A
dat_agg_diff <- 
  dat_diff %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75)
  ) %>%
  ungroup()
# combine
dat_out <- dat_agg %>% rbind(dat_agg_diff) 

# save 
saveRDS(dat_out, fpath_out)
rm(dat_out)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LMM: vary N_obs and error(s) sd 

fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-12-03-lmm_testcoef_vary_nobs_sd_raw")
fpath_out <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-12-03-lmm_testcoef_vary_nobs_sd_agg.rds")

# read and combine data 
fnames <- list.files(fdir_raw, full.names = TRUE)
fnames <- fnames[grepl("arrayjob", fnames)]
message(paste0("Raw files count: ", length(fnames)))
dat_l  <- lapply(fnames, readRDS)
dat    <- do.call("rbind", dat_l)
dim(dat)
head(dat)
table(dat$name)
table(dat$eff_tar, useNA = "ifany")
# aggregate data: estimated power  
dat_agg <- 
  dat %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma, sd_tau) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75),
    prop_success_mean = mean(prop_success),
    prop_success_min = min(prop_success),
    prop_success_max = max(prop_success)
  ) %>%
  ungroup()
summary(dat_agg$prop_success_mean)
summary(dat_agg$prop_success_min)
dat_agg <- dat_agg %>% select(-c(prop_success_mean, prop_success_min, prop_success_max))
dim(dat_agg)
table(dat_agg$eff_tar, useNA = "ifany")

# aggregate data: true mean 
dat_truepower_agg <- 
  dat %>%
  filter(name == "true_power") %>% 
  group_by(N_tar, N_obs, name, eff_tar, sd_sigma, sd_tau) %>%
  summarise(value = mean(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  ungroup() %>%
  select(-N_obs) %>%
  as.data.frame()
table(dat_truepower_agg$eff_tar, useNA = "ifany")
# compare with true power
dat_diff_A <- 
  dat %>%
  filter(name != "true_power") %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(dat_truepower_agg, by = c("N_tar", "eff_tar", "sd_sigma", "sd_tau")) %>%
  mutate(
    pe_ups_true =  100 * (upstrap_power - true_power)/true_power
  ) %>%
  select(-c(upstrap_power, true_power)) %>%
  pivot_longer(cols = starts_with("pe"))
dim(dat_diff_A)
dat_diff <- dat_diff_A
dat_agg_diff <- 
  dat_diff %>%
  group_by(N_tar, N_obs, name, eff_tru, eff_tar, sd_sigma, sd_tau) %>%
  summarise(
    cnt = sum(!is.na(value)),
    value_mean = mean(value),
    value_sd = sd(value),
    value_q25 = quantile(value, 0.25),
    value_q50 = median(value),
    value_q75 = quantile(value, 0.75)
  ) %>%
  ungroup()
# combine
dat_out <- dat_agg %>% rbind(dat_agg_diff) 

# save 
saveRDS(dat_out, fpath_out)
rm(dat_out)


# ------------------------------------------------------------------------------
#' cd $ups
#' git add --all
#' git commit -m 'add updated results from the cluster'
#' git push

