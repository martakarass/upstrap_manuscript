
rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
# source(here::here("numerical_experiments/R/config_figures.R"))
# source(here::here("numerical_experiments/R/config_utils.R"))

# dir to save results 
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-06-28-onesample_ttest_agg")
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-28-onesample_ttest_raw")
res_fdir_plt  <- paste0(here::here(), "/numerical_experiments/results_figures/2021-06-28-onesample_ttest")

# experiment parameters
N0 <- 45
N1_min <- 5
N1_max <- 200 
N1_grid <- N1_min : N1_max
N1_grid_l <- length(N1_grid)

# data generating model
mu0   <- 0.3
simga2 <- 1

# number of repetitions of experiment 
R_rep   <- 1000 
# number of boot repetitions within one experiment, one setup
B_boot  <- 1000

# effect of interest
mu1_grid <- c(0.2, 0.3, 0.4)


# ------------------------------------------------------------------------------
# READ DATA 

fnames_all   <- list.files(res_fdir_raw, full.names = TRUE)
dat <- do.call("rbind", lapply(fnames_all, readRDS))
dat <- dat %>% mutate(mu1 = ifelse(is.na(mu1), "as in observed sample x", mu1)) 
dim(dat)

head(dat)
tail(dat)


# ------------------------------------------------------------------------------
# ADD ASYMPTOTIC

for (mu1 in mu1_grid){
    value <- sapply(N1_grid, function(n_tmp){
        out_test <- power.t.test(n = n_tmp, delta = mu1, sd = sqrt(simga2), sig.level = 0.05, type = "one.sample", alternative = "two.sided")
        out_test$power 
    })
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(1, N1_grid_l)
    mat_out_tmp$name          <- "asymptotic_power"
    mat_out_tmp$mu0           <- mu0
    mat_out_tmp$mu1           <- mu1
    mat_out_tmp$value <- value
    dat <- rbind(dat, mat_out_tmp)
    rm(mat_out_tmp, value)
}
head(dat)
tail(dat)


# ------------------------------------------------------------------------------
# AGGREGATE

dat_agg <- 
    dat %>%
    group_by(N0, N1, mu0, mu1, name) %>%
    summarise(value_mean = mean(value),
              value_median = median(value)) %>%
    ungroup()

dat_agg_long <- 
    dat_agg %>% 
    pivot_longer(cols = starts_with("value_"), names_to = "agg_stat") %>%
    mutate(mu1_label = paste0("effect size = ", mu1))

dim(dat_agg_long)
head(dat_agg_long)
tail(dat_agg_long)
    

# ------------------------------------------------------------------------------
# PLOT

plt_ncol <- 2

plt1 <-
    ggplot(dat_agg_long %>% filter(name == "upstrap_power"), 
           aes(x = N1, y = value, color = agg_stat)) + 
    geom_hline(yintercept = 0.8, color = "blue", linetype = 2) + 
    geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
    geom_line() + 
    geom_line(data = dat_agg_long %>% filter(name == "asymptotic_power"),
              aes(x = N1, y = value), 
              color = "black") + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
    facet_wrap(~ mu1_label, ncol = plt_ncol) + 
    labs(x = "Target sample size M", y = "Power upstrap_i", color = "") + 
    theme(legend.position = c(0.80, 0.08))
plt1


plt2 <-
    ggplot(dat_agg_long %>% filter(name == "powerttest_power"), 
           aes(x = N1, y = value, color = agg_stat)) + 
    geom_hline(yintercept = 0.8, color = "blue", linetype = 2) + 
    geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
    geom_line() + 
    geom_line(data = dat_agg_long %>% filter(name == "asymptotic_power"),
              aes(x = N1, y = value), 
              color = "black") + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
    facet_wrap(~ mu1_label, ncol = plt_ncol) + 
    labs(x = "Target sample size M", y = "Power powerttest_i", color = "") + 
    theme(legend.position = c(0.80, 0.08))
plt2


plt3_df_a <- dat %>% filter(name == "upstrap_power") %>% rename(value_a = value) %>% select(-name)
plt3_df_b <- dat %>% filter(name == "powerttest_power") %>% rename(value_b = value) %>% select(-name)
plt3_df <- 
    plt3_df_a %>% 
    left_join(plt3_df_b, by = c("N1", "N0", "arrayjob_idx", "mu0", "mu1")) %>%
    mutate(value = value_a - value_b,
           name = "upstrap_powerttest_power_diff")
plt3_dat_agg_long <- 
    plt3_df %>%
    group_by(N0, N1, mu0, mu1, name) %>%
    summarise(value_mean = mean(value),
              value_median = median(value)) %>%
    ungroup() %>% 
    pivot_longer(cols = starts_with("value_"), names_to = "agg_stat") %>%
    mutate(mu1_label = paste0("effect size = ", mu1)) 

plt3 <-
    ggplot(plt3_dat_agg_long %>% filter(name == "upstrap_powerttest_power_diff"), 
           aes(x = N1, y = value, color = agg_stat)) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_vline(xintercept = N0, linetype = 2, size = 0.3) + 
    geom_line() + 
    # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + 
    facet_wrap(~ mu1_label, ncol = plt_ncol) + 
    labs(x = "Target sample size M", y = "Power diff (upstrap_i - powerttest_i)") + 
    theme(legend.position = c(0.80, 0.14))
plt3


plt_fpath <- paste0(res_fdir_plt, "/Power upstrap_i.png")
ggsave(filename = plt_fpath, plot = plt1, width = 6, height = 6)

plt_fpath <- paste0(res_fdir_plt, "/Power powerttest_i.png")
ggsave(filename = plt_fpath, plot = plt2, width = 6, height = 6)

plt_fpath <- paste0(res_fdir_plt, "/Power diff upstrap_i - powerttest_i.png")
ggsave(filename = plt_fpath, plot = plt3, width = 6, height = 6)







