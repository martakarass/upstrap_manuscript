
rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
# source(here::here("numerical_experiments/R/config_figures.R"))
# source(here::here("numerical_experiments/R/config_utils.R"))

# dir to save results 
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-06-28-lm_agg")
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-28-lm_raw")
res_fdir_plt  <- paste0(here::here(), "/numerical_experiments/results_figures/2021-06-28-lm")

N0 <- 31


# ------------------------------------------------------------------------------
# READ DATA 

fnames_all   <- list.files(res_fdir_raw, full.names = TRUE)
dat <- do.call("rbind", lapply(fnames_all, readRDS))
dat <- dat %>% mutate(effect1 = ifelse(is.na(effect1), "as in observed sample x", effect1)) 
dim(dat)

head(dat)
tail(dat)
table(dat$effect1, useNA = "ifany")
table(dat$N0, useNA = "ifany")
table(dat$N1, useNA = "ifany")
table(dat$name, useNA = "ifany")


# ------------------------------------------------------------------------------
# AGGREGATE

dat_agg <- 
    dat %>%
    group_by(N0, N1, effect0, effect1, name) %>%
    summarise(value_mean = mean(value),
              value_median = median(value)) %>%
    ungroup()

dat_agg_long <- 
    dat_agg %>% 
    pivot_longer(cols = starts_with("value_"), names_to = "agg_stat") %>%
    mutate(effect1_label = paste0("effect size = ", effect1))

dim(dat_agg_long)
head(dat_agg_long)
tail(dat_agg_long)
table(dat_agg_long$name)
    

# ------------------------------------------------------------------------------
# PLOT

plt1_df_asymp <- dat_agg_long %>% filter(name == "run_result", agg_stat == "value_mean")
plt1 <-
    ggplot(dat_agg_long %>% filter(name == "upstrap_power"), 
           aes(x = N1, y = value, color = agg_stat)) + 
    geom_hline(yintercept = 0.8, color = "blue", linetype = 2) + 
    geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
    geom_line() + 
    geom_line(data = plt1_df_asymp,
              aes(x = N1, y = value), 
              color = "black") + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
    facet_wrap(~ effect1_label, ncol = 2) + 
    labs(x = "Target sample size M", y = "Power upstrap_i", color = "") + 
    theme(legend.position = c(0.80, 0.08))
plt1


plt1_df_asymp <- dat_agg_long %>% filter(name == "run_result", agg_stat == "value_mean")
plt2 <-
    ggplot(dat_agg_long %>% filter(name == "simr_power"), 
           aes(x = N1, y = value, color = agg_stat)) + 
    geom_hline(yintercept = 0.8, color = "blue", linetype = 2) + 
    geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
    geom_line() + 
    geom_line(data = plt1_df_asymp,
              aes(x = N1, y = value), 
              color = "black") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
    facet_wrap(~ effect1_label, ncol = 2) + 
    labs(x = "Target sample size M", y = "Power powerttest_i", color = "") + 
    theme(legend.position = c(0.80, 0.08))
plt2


# plt3_df_a <- dat %>% filter(name == "upstrap_power") %>% rename(value_a = value) %>% select(-name)
# plt3_df_b <- dat %>% filter(name == "powerttest_power") %>% rename(value_b = value) %>% select(-name)
# plt3_df <- 
#     plt3_df_a %>% 
#     left_join(plt3_df_b, by = c("N1", "N0", "arrayjob_idx", "effect0", "effect1")) %>%
#     mutate(value = value_a - value_b,
#            name = "upstrap_powerttest_power_diff")
# plt3_dat_agg_long <- 
#     plt3_df %>%
#     group_by(N0, N1, effect0, effect1, name) %>%
#     summarise(value_mean = mean(value),
#               value_median = median(value)) %>%
#     ungroup() %>% 
#     pivot_longer(cols = starts_with("value_"), names_to = "agg_stat") %>%
#     mutate(effect1_label = paste0("effect size = ", effect1)) 
# 
# plt3 <-
#     ggplot(plt3_dat_agg_long %>% filter(name == "upstrap_powerttest_power_diff"), 
#            aes(x = N1, y = value, color = agg_stat)) + 
#     geom_hline(yintercept = 0, linetype = 2) + 
#     geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
#     geom_line() + 
#     # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
#     facet_wrap(~ effect1_label, ncol = 1) + 
#     labs(x = "Target sample size M", y = "Power diff (upstrap_i - powerttest_i)") + 
#     theme(legend.position = c(0.80, 0.14))
# plt3


plt_fpath <- paste0(res_fdir_plt, "/Power upstrap_i.png")
ggsave(filename = plt_fpath, plot = plt1, width = 6, height = 6)

plt_fpath <- paste0(res_fdir_plt, "/Power simrpower_i.png")
ggsave(filename = plt_fpath, plot = plt2,width = 6, height = 6)

# plt_fpath <- paste0(res_fdir_plt, "/Power diff upstrap_i - powerttest_i.png")
# ggsave(filename = plt_fpath, plot = plt3, width = 3.5, height = 8)







