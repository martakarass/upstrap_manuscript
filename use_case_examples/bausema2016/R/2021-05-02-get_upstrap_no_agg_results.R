#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# rm(list = ls())

#' For formatting the data, this script reuses code attached to the following paper: 
#' 
#' - Paper citation: Naudet F, Sakarovitch C, Janiaud P, Cristea I, Fanelli D, Moher D et al. Data sharing and reanalysis of randomized controlled trials in leading biomedical journals with a full data sharing policy: survey of studies published in The BMJ and PLOS Medicine BMJ 2018; 360 :k400 doi:10.1136/bmj.k400
#' - Paper URL: https://www.bmj.com/content/360/bmj.k400
#' - Code URL: https://osf.io/7ghfa/
#' 
#' Notes: 
#' cd $ups 
#' cd use_case_examples/bausema2016/R
#' Rnosave 2021-05-02-get_upstrap_no_agg_results.R -l mem_free=20G,h_vmem=20G,h_stack=256M -t 1-78 -tc 90 -N JOB_quality_flag
#' 
#' Rnosave 2021-05-02-get_upstrap_no_agg_results.R -t 1-78 -N JOB_bausema



# move the files needed 
# put /Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/use_case_examples/bausema2016/data_raw_csv/intervention_dataset_REDHOT_anonymous_updated.csv /users/mkaras/_PROJECTS/upstrap_manuscript/use_case_examples/bausema2016/data_raw_csv/intervention_dataset_REDHOT_anonymous_updated.csv


library(data.table)
library(dplyr)
library(matrixStats)
library(geepack)
options(scipen=999)

job_idx  <-  as.numeric(Sys.getenv("SGE_TASK_ID"))
# job_idx <- 10
set.seed(123)

# create dir directory if does not exist
fdir_tmp  <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-05-02-upstrap_no_agg_results/")
if (!dir.exists(fdir_tmp)) dir.create(fdir_tmp)


# ------------------------------------------------------------------------------
# READ AND FORMAT DATA

# read data with individual outcome status  
dat_fname <- paste0(here::here(), "/use_case_examples/bausema2016/data_raw_csv/intervention_dataset_REDHOT_anonymous_updated.csv")
dat0 <- fread(dat_fname) %>% as.data.frame()
dat0 <- dat0[!is.na(dat0$PCR_18S_result), ]
# add label for intervention/control cluster
cluster_arm <- c(
  "CONTROL","CONTROL","CONTROL","INTERVENTION","CONTROL",
  "INTERVENTION","INTERVENTION","INTERVENTION","CONTROL","INTERVENTION"
)
cluster_arm_df <- data.frame(cluster_arm = cluster_arm, cluster_number = 1:10)
dat <-  dat0 %>% left_join(cluster_arm_df, by = "cluster_number")
dat <- dat

# recode the data set according to Naudet's code 
# open eaves 
dat$eavetype[dat$eavetype=="Close"] <- "Closed"
dat$EAVE[dat$eavetype=="Open" | dat$eavetype=="Partially Open" ] <- "Partially Open"
dat$EAVE[dat$eavetype=="Closed" ] <- "Closed"
dat$EAVE[dat$eavetype=="" ] <- NA
# age group
dat$AGEGROUP[as.numeric(as.character(dat$age))<6] <- "<=5"
dat$AGEGROUP[as.numeric(as.character(dat$age))>=6 & as.numeric(as.character(dat$age))<16] <- "6-15"
dat$AGEGROUP[as.numeric(as.character(dat$age))>=16 & as.numeric(as.character(dat$age))<26] <- "16-25"
dat$AGEGROUP[as.numeric(as.character(dat$age))>=26] <- ">25"
# altitude group
dat$ALTITUDE[dat$elevation_group=="1350-1449" | dat$elevation_group=="1396-1424" | dat$elevation_group=="1425-1449"] <- "INF1450"
dat$ALTITUDE[dat$elevation_group=="1450-1474" | dat$elevation_group=="1450-1499" | dat$elevation_group=="1475-1499"] <- "1450-1500"
dat$ALTITUDE[dat$elevation_group=="1550-1605" | dat$elevation_group=="1500+" | dat$elevation_group=="1500-1549"] <- "SUP1500"

# rename for convenience
dat$y <- dat$PCR_18S_result
dat$zone <- as.character(dat$distance_hotspot_categories)
dat$zone[dat$zone == "0"] <- "hotspot"
dat$zone[dat$zone == "1"] <- "eval_zone_1"
dat$zone[dat$zone == "2"] <- "eval_zone_2"
dat$time_point <- as.character(dat$survey_round)
dat$time_point[dat$time_point == "1"] <- "baseline"
dat$time_point[dat$time_point == "2"] <- "wk_8"
dat$time_point[dat$time_point == "3"] <- "wk_16"


# ------------------------------------------------------------------------------
# COMPUTE ADJUSTED OUTCOME (RESIDUALS)

# get the prevalence baseline -- cluster- and zone-specific  
dat_agg_baseline <- dat %>% 
  filter(time_point == "baseline") %>%
  group_by(cluster_number, zone) %>%
  summarise(y = mean(y, na.rm = TRUE)) %>%
  rename(y_BASELINE = y) %>%
  as.data.frame()
dat <- dat %>% left_join(dat_agg_baseline, by = c("zone", "cluster_number"))
# make empty vector of predicted values
dat$y_PRED <- NA

# do logistic regression to predict the outcome -- response time-specific ------
# - without using intervention arm 
# - we model per time_point
# - use type = "response" that gives the predicted probabilities. 
# time_point == wk_8
dat_rowidx <- which(dat$time_point == "wk_8")
dat_tmp <- dat[dat_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$y) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$y_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
dat$y_PRED[dat_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
# mean(dat$y_PRED[dat_rowidx], na.rm = TRUE)
# [1] 0.1316203

# time_point == wk_16
dat_rowidx <- which(dat$time_point == "wk_16")
dat_tmp <- dat[dat_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$y) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$y_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
dat$y_PRED[dat_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
# mean(dat$y_PRED[dat_rowidx], na.rm = TRUE)
# [1] 0.1122656

# generate final data frame used in fits 
fit_df_all <- 
  dat %>%
  filter(!is.na(y_PRED)) %>%
  mutate(cluster_arm = factor(as.character(cluster_arm), levels = c("CONTROL", "INTERVENTION")))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# single fits

fit_GEE_formula <- formula(y ~ y_PRED + cluster_arm)

# wk_8, hotspot  ---------------------------------------------------------------
dat_sub <- fit_df_all %>% filter(time_point == "wk_8", zone == "hotspot")
fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                  data = dat_sub, id = dat_sub$cluster_number, corstr = "exchangeable")
fit_GEE_s <- summary(fit_GEE)
# fit_GEE_s
# round(fit_GEE_s$coefficients$Estimate, 3)
# round(exp(fit_GEE_s$coefficients$Estimate), 3)
# round(fit_GEE_s$coefficients$Std.err, 3)
# round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# wk_8, eval_zone_1  ---------------------------------------------------------------
dat_sub <- fit_df_all %>% filter(time_point == "wk_8", zone == "eval_zone_1")
fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                  data = dat_sub, id = dat_sub$cluster_number, corstr = "exchangeable")
fit_GEE_s <- summary(fit_GEE)
# fit_GEE_s
# round(fit_GEE_s$coefficients$Estimate, 3)
# round(exp(fit_GEE_s$coefficients$Estimate), 3)
# round(fit_GEE_s$coefficients$Std.err, 3)
# round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# wk_8, eval_zone_2  ---------------------------------------------------------------
dat_sub <- fit_df_all %>% filter(time_point == "wk_8", zone == "eval_zone_2")
fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                  data = dat_sub, id = dat_sub$cluster_number, corstr = "exchangeable")
fit_GEE_s <- summary(fit_GEE)
# fit_GEE_s
# round(fit_GEE_s$coefficients$Estimate, 3)
# round(exp(fit_GEE_s$coefficients$Estimate), 3)
# round(fit_GEE_s$coefficients$Std.err, 3)
# round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# ------------------------------------------------------------------------------
# DO UPSTRAP

M_min    <- 5
M_max    <- 30
M_grid   <- M_min : M_max
M_grid_l <- length(M_grid)
B_boot   <- 1000 * 10

# define grid of upstrap analysis parameters
param_df <- expand.grid(
  time_point = c("wk_8"),
  zone = c("hotspot", "eval_zone_1", "eval_zone_2"),
  sample_size_M = M_grid
)
# dim(param_df)
# [1] 78  3

# pull upstrap analysis parameters specific to current job 
zone_tmp       <- param_df[job_idx, "zone"]
time_point_tmp <- param_df[job_idx, "time_point"]
M_tmp          <- param_df[job_idx, "sample_size_M"]
message(paste0("time_point: ", time_point_tmp, ", zone_tmp: ", zone_tmp, ", M_tmp: ", M_tmp))

# generate vectors with cluster numbers for two arms 
cluster_meta <- fit_df_all %>% select(cluster_number, cluster_arm) %>% distinct()
cn_C <- cluster_meta %>% filter(cluster_arm == "CONTROL") %>% pull(cluster_number)
cn_I <- cluster_meta %>% filter(cluster_arm == "INTERVENTION") %>% pull(cluster_number)

t1 <- Sys.time()
# get data subset specific to current job  
dat_sub <- fit_df_all %>% filter(time_point == time_point_tmp, zone == zone_tmp)

# sample row index for current upstrap sampling: control, intervention 
boot_cn_C <- matrix(sample(x = cn_C, size = (B_boot * M_tmp), replace = TRUE), nrow = B_boot)
boot_cn_I <- matrix(sample(x = cn_I, size = (B_boot * M_tmp), replace = TRUE), nrow = B_boot)

# model formula
fit_GEE_formula <- formula(y ~ y_PRED + cluster_arm)

# repeat procedure B_boot-times 
boot_power_vec <- numeric(B_boot)
message("Starting the upstrap procedure...")
t1 <- Sys.time()
for (j in 1 : B_boot){ # j <- 1
  tryCatch({
    boot_idx_C <- boot_cn_C[j, ]
    boot_idx_I <- boot_cn_I[j, ]
    boot_idx   <- c(boot_idx_C, boot_idx_I) 
    # make a data set that corresponds to current boot resampling 
    boot_data_list <- lapply(1 : length(boot_idx), function(el) {
      dat_sub %>% filter(cluster_number == boot_idx[el]) %>%
        mutate(cluster_number = el)
    })
    boot_data <- do.call("rbind", boot_data_list)
    # fit the model 
    fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                      data = boot_data, id = boot_data$cluster_number, corstr = "exchangeable")
    fit_GEE_s <- summary(fit_GEE)
    # store information 
    boot_power_vec[j] <- fit_GEE_s$coefficients[3, 4] < 0.05 
  }, error = function(e) {
    message(e)
  })
  # print time message
  if (j %% 100 == 0) {
    message(paste0("j: ", j, " Time elapsed: ", round(as.numeric(Sys.time() - t1, units = "hours"), 3), " hours"))
    # append results to data frame
    out_df <- data.frame(
      time_point = time_point_tmp, 
      zone = zone_tmp, 
      sample_size_M = M_tmp,
      power = mean(boot_power_vec, na.rm = TRUE),
      B_boot = B_boot,
      B_boot_success = sum(!is.na(boot_power_vec)),
      eval_time_secs = round(as.numeric(Sys.time() - t1, units = "secs"))
    )
    # save results to file 
    fpath_tmp <- paste0(fdir_tmp, job_idx, ".rds")
    saveRDS(out_df, fpath_tmp)
  }
}
message("Upstrap procedure completed.")

# append results to data frame
out_df <- data.frame(
  time_point = time_point_tmp, 
  zone = zone_tmp, 
  sample_size_M = M_tmp,
  power = mean(boot_power_vec, na.rm = TRUE),
  B_boot = B_boot,
  B_boot_success = sum(!is.na(boot_power_vec)),
  eval_time_secs = round(as.numeric(Sys.time() - t1, units = "secs"))
)

# save results to file 
fpath_tmp <- paste0(fdir_tmp, job_idx, ".rds")
saveRDS(out_df, fpath_tmp)
message(paste0("Results saved. Time elapsed: ", round(as.numeric(Sys.time() - t1, units = "hours"), 3), " hours"))
