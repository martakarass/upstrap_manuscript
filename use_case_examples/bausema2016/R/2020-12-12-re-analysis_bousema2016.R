
rm(list = ls())
library(data.table)
library(dplyr)

dat_fname <- paste0(here::here(), "/use_case_examples/data_processed/bousema2016/intervention_dataset_REDHOT_anonymous_updated.csv")
dat0 <- fread(dat_fname) %>% as.data.frame()
names(dat0)
dim(dat0)
# [1] 12989    48
sum(!is.na(dat0$PCR_18S_result))
# [1] 12938
table(dat0$PCR_18S_result)
dat1 <- dat0[!is.na(dat0$PCR_18S_result), ]
dim(dat1)
# [1] 12938    48

# definition of intervention/control cluster
cluster_arm <- c("CONTROL","CONTROL","CONTROL","INTERVENTION","CONTROL","INTERVENTION","INTERVENTION","INTERVENTION","CONTROL","INTERVENTION")
cluster_arm_df <- data.frame(cluster_arm = cluster_arm, cluster_number = 1:10)
dat2 <-  dat1 %>% left_join(cluster_arm_df, by = "cluster_number")
  

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# #1 -- REPLICATE UNDJUSTED ANALYSIS, AS IN THE PAPER 

TAB1 <- dat2 
summary(as.factor(TAB1$PCR_18S_result))
summary(as.factor(TAB1$cluster_number))
summary(as.factor(TAB1$survey_round))
summary(as.factor(TAB1$hotspot_evaluation))
summary(as.factor(TAB1$distance_hotspot_categories))
summary(as.factor(TAB1$cluster_arm))
table(TAB1$hotspot_evaluation, TAB1$distance_hotspot_categories)
table(TAB1$cluster_number, TAB1$cluster_arm)

# aggregate to get proportion of cases in each (cluster, cluster_zone)
TAB1_agg <- TAB1 %>% 
  group_by(cluster_number, distance_hotspot_categories, cluster_arm, survey_round) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  as.data.frame()
dim(TAB1_agg)
# [1] 90  5

# distance_hotspot_categories == 0 (hotspot) -----------------------------------
TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 0, survey_round == 2)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 0, survey_round == 3)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

# distance_hotspot_categories == 1 (evaluation zone, closer) -----------------------------------
TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 1, survey_round == 2)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 1, survey_round == 3)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

# distance_hotspot_categories == 2 (evaluation zone, further) -----------------------------------
TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 2, survey_round == 2)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

TAB1_agg_sub <- TAB1_agg %>% filter(distance_hotspot_categories == 2, survey_round == 3)
t.test(
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result"],
  TAB1_agg_sub[TAB1_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result"]
)

rm(TAB1, TAB1_agg, TAB1_agg_sub)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# #2 -- REPLICATE UNDJUSTED ANALYSIS, BUT USE DIFF WITH BASELINE 

# get diff baseline 
TAB2 <- dat2
TAB2_agg0 <- TAB2 %>% 
  group_by(cluster_number, distance_hotspot_categories, cluster_arm, survey_round) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  as.data.frame()
TAB2_agg_baseline <- TAB2 %>% 
  filter(survey_round == 1) %>%
  group_by(cluster_number, distance_hotspot_categories) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  rename(PCR_18S_result_BASELINE = PCR_18S_result) %>%
  as.data.frame()
TAB2_agg <- TAB2_agg0 %>% 
  left_join(TAB2_agg_baseline, by = c("cluster_number", "distance_hotspot_categories")) %>%
  mutate(PCR_18S_result_DIFF = PCR_18S_result - PCR_18S_result_BASELINE)


# distance_hotspot_categories == 0 (hotspot) -----------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 0, survey_round == 2)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 0, survey_round == 3)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

# distance_hotspot_categories == 1 (evaluation zone, closer) -----------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 1, survey_round == 2)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 1, survey_round == 3)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

# distance_hotspot_categories == 2 (evaluation zone, further) -----------------------------------
TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 2, survey_round == 2)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

TAB2_agg_sub <- TAB2_agg %>% filter(distance_hotspot_categories == 2, survey_round == 3)
t.test(
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_DIFF"],
  TAB2_agg_sub[TAB2_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_DIFF"]
)

rm(TAB2, TAB2_agg0, TAB2_agg, TAB2_agg_sub, TAB2_agg_baseline)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# #3 -- REPLICATE ADJUSTED ANALYSIS, AS IN THE PAPER

rm(TAB3)
TAB3 <- dat2 

# recode the data set according to Naudet's code -------------------------------
# open eaves 
TAB3$eavetype[TAB3$eavetype=="Close"] <- "Closed"
TAB3$EAVE[TAB3$eavetype=="Open" | TAB3$eavetype=="Partially Open" ] <- "Partially Open"
TAB3$EAVE[TAB3$eavetype=="Closed" ] <- "Closed"
TAB3$EAVE[TAB3$eavetype=="" ] <- NA
# age group
TAB3$AGEGROUP[as.numeric(as.character(TAB3$age))<6] <- "<=5"
TAB3$AGEGROUP[as.numeric(as.character(TAB3$age))>=6 & as.numeric(as.character(TAB3$age))<16] <- "6-15"
TAB3$AGEGROUP[as.numeric(as.character(TAB3$age))>=16 & as.numeric(as.character(TAB3$age))<26] <- "16-25"
TAB3$AGEGROUP[as.numeric(as.character(TAB3$age))>=26] <- ">25"
# altitude group
TAB3$ALTITUDE[TAB3$elevation_group=="1350-1449" | TAB3$elevation_group=="1396-1424" | TAB3$elevation_group=="1425-1449"] <- "INF1450"
TAB3$ALTITUDE[TAB3$elevation_group=="1450-1474" | TAB3$elevation_group=="1450-1499" | TAB3$elevation_group=="1475-1499"] <- "1450-1500"
TAB3$ALTITUDE[TAB3$elevation_group=="1550-1605" | TAB3$elevation_group=="1500+" | TAB3$elevation_group=="1500-1549"] <- "SUP1500"

summary(as.factor(TAB3$EAVE))
summary(as.factor(TAB3$AGEGROUP))
summary(as.factor(TAB3$ALTITUDE))


# get the prevalence baseline -- cluster- and zone-specific  -------------------
TAB3_agg_baseline <- TAB3 %>% 
  filter(survey_round == 1) %>%
  group_by(cluster_number, distance_hotspot_categories) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  rename(PCR_18S_result_BASELINE = PCR_18S_result) %>%
  as.data.frame()
TAB3 <- TAB3 %>% left_join(TAB3_agg_baseline, by = c("distance_hotspot_categories", "cluster_number"))
# make empty vector of predicted values
TAB3$PCR_18S_result_PRED <- NA


# do logistic regression to predict the outcome -- response time-specific ------
# - without using intervention arm 
# - we model per survey_round
# - use type = "response" that gives the predicted probabilities. 

# survey_round == 2
TAB3_rowidx <- which(TAB3$survey_round == 2)
dat_tmp <- TAB3[TAB3_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$PCR_18S_result) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$PCR_18S_result_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
TAB3$PCR_18S_result_PRED[TAB3_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
mean(TAB3$PCR_18S_result_PRED[TAB3_rowidx], na.rm = TRUE)
# [1] 0.1316203

# survey_round == 3
TAB3_rowidx <- which(TAB3$survey_round == 3)
dat_tmp <- TAB3[TAB3_rowidx, ]
MODEL_tmp <- glm(
  as.factor(dat_tmp$PCR_18S_result) ~ dat_tmp$AGEGROUP + as.factor(dat_tmp$gender) + dat_tmp$ALTITUDE + dat_tmp$EAVE + dat_tmp$PCR_18S_result_BASELINE, 
  binomial(link = 'logit'), na.action = na.exclude)
TAB3$PCR_18S_result_PRED[TAB3_rowidx] <-  predict(MODEL_tmp, type = "response", na.action = na.exclude)
mean(TAB3$PCR_18S_result_PRED[TAB3_rowidx], na.rm = TRUE)
# [1] 0.1122656


# compute the residuals between predicted and actual outcomes ------------------
# (act-exp) which matches the convention at wiki 
mean(TAB3$PCR_18S_result_PRED, na.rm = TRUE)
hist(TAB3$PCR_18S_result_PRED, breaks = 50)

TAB3_agg <- TAB3 %>% 
  group_by(cluster_number, distance_hotspot_categories, cluster_arm, survey_round) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE),
            PCR_18S_result_PRED = mean(PCR_18S_result_PRED, na.rm = TRUE)) %>%
  mutate(PCR_18S_result_RESID = PCR_18S_result - PCR_18S_result_PRED) %>%
  as.data.frame()
TAB3_agg

# distance_hotspot_categories == 0 (hotspot) -----------------------------------
TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 0, survey_round == 2)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 0, survey_round == 3)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

# distance_hotspot_categories == 1 (evaluation zone, closer) -----------------------------------
TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 1, survey_round == 2)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 1, survey_round == 3)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

# distance_hotspot_categories == 2 (evaluation zone, further) -----------------------------------
TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 2, survey_round == 2)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

TAB3_agg_sub <- TAB3_agg %>% filter(distance_hotspot_categories == 2, survey_round == 3)
t.test(
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "CONTROL", "PCR_18S_result_RESID"],
  TAB3_agg_sub[TAB3_agg_sub$cluster_arm == "INTERVENTION", "PCR_18S_result_RESID"]
)

