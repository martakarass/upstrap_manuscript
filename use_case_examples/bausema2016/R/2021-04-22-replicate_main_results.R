
rm(list = ls())
library(data.table)
library(dplyr)

#' This script reuses code attached to the following paper: 
#' 
#' - Paper citation: Naudet F, Sakarovitch C, Janiaud P, Cristea I, Fanelli D, Moher D et al. Data sharing and reanalysis of randomized controlled trials in leading biomedical journals with a full data sharing policy: survey of studies published in The BMJ and PLOS Medicine BMJ 2018; 360 :k400 doi:10.1136/bmj.k400
#' - Paper URL: https://www.bmj.com/content/360/bmj.k400
#' - Code URL: https://osf.io/7ghfa/

# make placeholder for paper's table
# df_tmp <- matrix(NA, nrow = (3 + 2 + 5 + 5), ncol = 5) %>% as.data.frame()
# stargazer::stargazer(df_tmp, summary = FALSE)


# ------------------------------------------------------------------------------
# read data 

# read data with individual outcome status  
dat_fname <- paste0(here::here(), "/use_case_examples/bausema2016/data_raw_csv/intervention_dataset_REDHOT_anonymous_updated.csv")
dat0 <- fread(dat_fname) %>% as.data.frame()
dat0 <- dat0[!is.na(dat0$PCR_18S_result), ]

table(dat0$PCR_18S_result)
dim(dat0)

# add label for intervention/control cluster
cluster_arm <- c(
  "CONTROL","CONTROL","CONTROL","INTERVENTION","CONTROL",
  "INTERVENTION","INTERVENTION","INTERVENTION","CONTROL","INTERVENTION"
)
cluster_arm_df <- data.frame(cluster_arm = cluster_arm, cluster_number = 1:10)
dat <-  dat0 %>% left_join(cluster_arm_df, by = "cluster_number")
  


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PART 1: UNADJSTED ANALYSIS 

TAB1 <- dat 
summary(as.factor(TAB1$PCR_18S_result))
summary(as.factor(TAB1$cluster_number))
summary(as.factor(TAB1$survey_round))
summary(as.factor(TAB1$hotspot_evaluation))
summary(as.factor(TAB1$distance_hotspot_categories))
summary(as.factor(TAB1$cluster_arm))
table(TAB1$hotspot_evaluation, TAB1$distance_hotspot_categories)
table(TAB1$cluster_number, TAB1$cluster_arm)

# rename for convenience
TAB1$y <- TAB1$PCR_18S_result
TAB1$zone <- as.character(TAB1$distance_hotspot_categories)
TAB1$zone[TAB1$zone == "0"] <- "hotspot"
TAB1$zone[TAB1$zone == "1"] <- "eval_zone_1"
TAB1$zone[TAB1$zone == "2"] <- "eval_zone_2"
TAB1$time_point <- as.character(TAB1$survey_round)
TAB1$time_point[TAB1$time_point == "1"] <- "baseline"
TAB1$time_point[TAB1$time_point == "2"] <- "wk_8"
TAB1$time_point[TAB1$time_point == "3"] <- "wk_16"

# aggregate to get proportion of cases in each (cluster, cluster_zone)
TAB1_agg <- TAB1 %>% 
  group_by(cluster_number, cluster_arm, zone, time_point) %>%
  summarise(y = mean(y, na.rm = TRUE)) %>%
  as.data.frame()
TAB1_agg


# ------------------------------------------------------------------------------
# baseline, hotspot  -----------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "baseline", zone == "hotspot")
TAB1_agg_sub %>% group_by(cluster_arm) %>% summarize(y = mean(y)) %>% mutate(y = round(y * 100, 1))

# baseline, eval_zone_1  -------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "baseline", zone == "eval_zone_1")
TAB1_agg_sub %>% group_by(cluster_arm) %>% summarize(y = mean(y)) %>% mutate(y = round(y * 100, 1))

# baseline, eval_zone_2  -------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "baseline", zone == "eval_zone_2")
TAB1_agg_sub %>% group_by(cluster_arm) %>% summarize(y = mean(y)) %>% mutate(y = round(y * 100, 1))


# ------------------------------------------------------------------------------
# week 8, hotspot  -------------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "wk_8", zone == "hotspot")
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "y"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "y"]
)

# week 8, eval_zone_1  -------------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "wk_8", zone == "eval_zone_1")
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "y"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "y"]
)

# week 8, eval_zone_2  -------------------------------------------------------------

TAB1_agg_sub <- TAB1_agg %>% filter(time_point == "wk_8", zone == "eval_zone_2")
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "y"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "y"]
)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PART 2: ADJUSTED ANALYSIS 

rm(TAB1)
TAB2 <- dat

# recode the data set according to Naudet's code -------------------------------
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

summary(as.factor(TAB2$EAVE))
summary(as.factor(TAB2$AGEGROUP))
summary(as.factor(TAB2$ALTITUDE))

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


# get the prevalence baseline -- cluster- and zone-specific  -------------------
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


# compute the residuals between predicted and actual outcomes ------------------
# (act-exp) which matches the convention at wiki 
TAB2_agg <- TAB2 %>% 
  group_by(cluster_number, cluster_arm, zone, time_point) %>%
  summarise(y = mean(y, na.rm = TRUE),
            y_PRED = mean(y_PRED, na.rm = TRUE)) %>%
  mutate(y_RESID = y - y_PRED) %>%
  as.data.frame()


# wk_8, hotspot  ---------------------------------------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(time_point == "wk_8", zone == "hotspot")
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "y_RESID"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "y_RESID"]
)

# wk_8, eval_zone_1  ---------------------------------------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(time_point == "wk_8", zone == "eval_zone_1")
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "y_RESID"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "y_RESID"]
)

# wk_8, eval_zone_2  ---------------------------------------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(time_point == "wk_8", zone == "eval_zone_2")
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "y_RESID"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "y_RESID"]
)


