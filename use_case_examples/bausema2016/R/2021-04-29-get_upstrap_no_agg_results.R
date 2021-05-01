
#' For formatting the data, this script reuses code attached to the following paper: 
#' 
#' - Paper citation: Naudet F, Sakarovitch C, Janiaud P, Cristea I, Fanelli D, Moher D et al. Data sharing and reanalysis of randomized controlled trials in leading biomedical journals with a full data sharing policy: survey of studies published in The BMJ and PLOS Medicine BMJ 2018; 360 :k400 doi:10.1136/bmj.k400
#' - Paper URL: https://www.bmj.com/content/360/bmj.k400
#' - Code URL: https://osf.io/7ghfa/

rm(list = ls())
library(data.table)
library(dplyr)
library(matrixStats)
library(geepack)
options(scipen=999)

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
TAB2 <- dat

# recode the data set according to Naudet's code 
# open eaves 
TAB2$eavetype[TAB2$eavetype=="Close"] <- "Closed"
TAB2$EAVE[TAB2$eavetype=="Open" | TAB2$eavetype=="Partially Open" ] <- "Partially Open"
TAB2$EAVE[TAB2$eavetype=="Closed" ] <- "Closed"
TAB2$EAVE[TAB2$eavetype=="" ] <- NA
# age group
TAB2$AGEGROUP[as.numeric(as.character(TAB2$age))<6] <- "<=5"
TAB2$AGEGROUP[as.numeric(as.character(TAB2$age))>=6 & as.numeric(as.character(TAB2$age))<16] <- "6-15"
TAB2$AGEGROUP[as.numeric(as.character(TAB2$age))>=16 & as.numeric(as.character(TAB2$age))<26] <- "16-25"
TAB2$AGEGROUP[as.numeric(as.character(TAB2$age))>=26] <- ">25"
# altitude group
TAB2$ALTITUDE[TAB2$elevation_group=="1350-1449" | TAB2$elevation_group=="1396-1424" | TAB2$elevation_group=="1425-1449"] <- "INF1450"
TAB2$ALTITUDE[TAB2$elevation_group=="1450-1474" | TAB2$elevation_group=="1450-1499" | TAB2$elevation_group=="1475-1499"] <- "1450-1500"
TAB2$ALTITUDE[TAB2$elevation_group=="1550-1605" | TAB2$elevation_group=="1500+" | TAB2$elevation_group=="1500-1549"] <- "SUP1500"

# rename for convenience
TAB2$y <- TAB2$PCR_18S_result
TAB2$zone <- as.character(TAB2$distance_hotspot_categories)
TAB2$zone[TAB2$zone == "0"] <- "hotspot"
TAB2$zone[TAB2$zone == "1"] <- "eval_zone_1"
TAB2$zone[TAB2$zone == "2"] <- "eval_zone_2"
TAB2$time_point <- as.character(TAB2$survey_round)
TAB2$time_point[TAB2$time_point == "1"] <- "baseline"
TAB2$time_point[TAB2$time_point == "2"] <- "wk_8"
TAB2$time_point[TAB2$time_point == "3"] <- "wk_16"


# ------------------------------------------------------------------------------
# COMPUTE ADJUSTED OUTCOME (RESIDUALS)

# get the prevalence baseline -- cluster- and zone-specific  
TAB2_agg_baseline <- TAB2 %>% 
  filter(time_point == "baseline") %>%
  group_by(cluster_number, zone) %>%
  summarise(y = mean(y, na.rm = TRUE)) %>%
  rename(y_BASELINE = y) %>%
  as.data.frame()
TAB2 <- TAB2 %>% left_join(TAB2_agg_baseline, by = c("zone", "cluster_number"))
# make empty vector of predicted values
TAB2$y_PRED <- NA

# do logistic regression to predict the outcome -- response time-specific ------
# - without using intervention arm 
# - we model per time_point
# - use type = "response" that gives the predicted probabilities. 
# time_point == wk_8
TAB2_rowidx <- which(TAB2$time_point == "wk_8")
dat_tmp <- TAB2[TAB2_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$y) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$y_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
TAB2$y_PRED[TAB2_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
mean(TAB2$y_PRED[TAB2_rowidx], na.rm = TRUE)
# [1] 0.1316203

# time_point == wk_16
TAB2_rowidx <- which(TAB2$time_point == "wk_16")
dat_tmp <- TAB2[TAB2_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$y) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$y_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
TAB2$y_PRED[TAB2_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
mean(TAB2$y_PRED[TAB2_rowidx], na.rm = TRUE)
# [1] 0.1122656


# investigate the aggregates  --------------------------------------------------

TAB2 %>% 
  group_by(time_point, zone) %>%
  summarise(
    y_PRED_NA_cnt = sum(is.na(y_PRED)),
    y_PRED_nonNA_cnt = sum(!is.na(y_PRED))
  )

TAB2 %>% 
  filter(!is.na(y_PRED)) %>%
  group_by(time_point, zone, cluster_arm) %>%
  summarise(
    y_mean = mean(y),
    y_PRED_mean = mean(y_PRED)
  )

fit_df_all <- 
  TAB2 %>%
  filter(!is.na(y_PRED)) %>%
  mutate(cluster_arm = factor(as.character(cluster_arm), levels = c("CONTROL", "INTERVENTION")))

fit_df_all %>% 
  group_by(time_point, zone, cluster_number, cluster_arm) %>%
  summarise(cnt = n())  %>%
  as.data.frame()


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
fit_GEE_s
round(fit_GEE_s$coefficients$Estimate, 3)
round(exp(fit_GEE_s$coefficients$Estimate), 3)
round(fit_GEE_s$coefficients$Std.err, 3)
round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# wk_8, eval_zone_1  ---------------------------------------------------------------
dat_sub <- fit_df_all %>% filter(time_point == "wk_8", zone == "eval_zone_1")
fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                  data = dat_sub, id = dat_sub$cluster_number, corstr = "exchangeable")
fit_GEE_s <- summary(fit_GEE)
fit_GEE_s
round(fit_GEE_s$coefficients$Estimate, 3)
round(exp(fit_GEE_s$coefficients$Estimate), 3)
round(fit_GEE_s$coefficients$Std.err, 3)
round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# wk_8, eval_zone_2  ---------------------------------------------------------------
dat_sub <- fit_df_all %>% filter(time_point == "wk_8", zone == "eval_zone_2")
fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                  data = dat_sub, id = dat_sub$cluster_number, corstr = "exchangeable")
fit_GEE_s <- summary(fit_GEE)
fit_GEE_s
round(fit_GEE_s$coefficients$Estimate, 3)
round(exp(fit_GEE_s$coefficients$Estimate), 3)
round(fit_GEE_s$coefficients$Std.err, 3)
round(fit_GEE_s$coefficients$`Pr(>|W|)`, 5)


# ------------------------------------------------------------------------------
# DO UPSTRAP

M_min    <- 5
M_max    <- 30
M_grid   <- M_min : M_max
M_grid_l <- length(M_grid)
B_boot   <- 10 * 10

# define objects to store results
time_point_out <- numeric()
zone_out       <- numeric()
power_out      <- numeric()
sample_size_M_out <- numeric()

cluster_meta <- fit_df_all %>% select(cluster_number, cluster_arm) %>% distinct()
cn_C <- cluster_meta %>% filter(cluster_arm == "CONTROL") %>% pull(cluster_number)
cn_I <- cluster_meta %>% filter(cluster_arm == "INTERVENTION") %>% pull(cluster_number)

t1 <- Sys.time()
set.seed(123)
# iterate over statistical tests (data subsets)
for (time_point_tmp in c("wk_8")){
  for (zone_tmp in c("hotspot", "eval_zone_1", "eval_zone_2")){ 
    # time_point_tmp <- "wk_8"; zone_tmp <- "hotspot"
    dat_sub <- fit_df_all %>% filter(time_point == time_point_tmp, zone == zone_tmp)
    
    # iterate over sample size M 
    for (M_tmp in rev(M_grid)){ # M_tmp <- 20
      message(paste0("time_point: ", time_point_tmp, ", zone_tmp: ", zone_tmp, ", M: ", M_tmp))
      # sample row index for current upstrap sampling: control, intervention 
      boot_cn_C <- matrix(sample(x = cn_C, size = (B_boot * M_tmp), replace = TRUE), nrow = B_boot)
      boot_cn_I <- matrix(sample(x = cn_I, size = (B_boot * M_tmp), replace = TRUE), nrow = B_boot)
      # repeat procedure B_boot-times 
      boot_power_vec <- numeric(B_boot)
      for (j in 1 : B_boot){ # j <- 1
        if (j %% 1000 == 0) print(j)
        boot_idx_C <- boot_cn_C[j, ]
        boot_idx_I <- boot_cn_I[j, ]
        boot_idx  <- c(boot_idx_C, boot_idx_I) 
        # make a data set that corresponds to 
        boot_data_list <- lapply(1 : length(boot_idx), function(el) {
          dat_sub %>% filter(cluster_number == boot_idx[el]) %>%
            mutate(cluster_number = el)
        })
        boot_data <- do.call("rbind", boot_data_list)
        fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                          data = boot_data, id = boot_data$cluster_number, corstr = "exchangeable")
        fit_GEE_s <- summary(fit_GEE)
        boot_power_vec[j] <- fit_GEE_s$coefficients[3, 4] < 0.05 
      }
      # store M-specific results
      time_point_out <- c(time_point_out, time_point_tmp)
      zone_out       <- c(zone_out, zone_tmp)
      power_out      <- c(power_out, mean(boot_power_vec, na.rm = TRUE))
      sample_size_M_out <- c(sample_size_M_out, M_tmp)
      t_tmp <- Sys.time()
      message(paste0("zone_tmp: ", zone_tmp, ", M_tmp: ", M_tmp, ", time elapsed: ", t_tmp - t1))
    }
    
  }
}
t2 <- Sys.time()
t2 - t1

# append results to data frame
out_df <- data.frame(
  time_point = time_point_out, 
  zone = zone_out, 
  power = power_out,
  sample_size_M = sample_size_M_out)
out_df

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/bausema2016/results/2021-04-30-upstrap_no_agg_results.rds")
saveRDS(object = out_df, file = out_power_df_fpath)








