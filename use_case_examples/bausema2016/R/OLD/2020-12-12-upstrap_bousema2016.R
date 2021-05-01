
#' This script upstraps data to estimate power

rm(list = ls())
library(data.table)
library(dplyr)
library(matrixStats)

# definition of intervention/control cluster
cluster_arm <- c("CONTROL","CONTROL","CONTROL","INTERVENTION","CONTROL",
                 "INTERVENTION","INTERVENTION","INTERVENTION","CONTROL","INTERVENTION")
cluster_arm_df <- data.frame(cluster_arm = cluster_arm, cluster_number = 1:10)

dat_fname <- paste0(here::here(), "/use_case_examples/data_processed/bousema2016/intervention_dataset_REDHOT_anonymous_updated.csv")
dat <- fread(dat_fname) %>% 
  as.data.frame() %>% 
  filter(!is.na(PCR_18S_result)) %>% 
  left_join(cluster_arm_df, by = "cluster_number")
dim(dat)
# [1] 12938    49


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PREPARE DATA SETS

# #1 -- DATA SET FOR  UNDJUSTED ANALYSIS (AS IN THE PAPER) ---------------------

TAB1 <- dat 
TAB1_agg <- TAB1 %>% 
  group_by(cluster_number, distance_hotspot_categories, cluster_arm, survey_round) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  as.data.frame()


# #2 -- DATA SET FOR ADJUSTED ANALYSIS (AS IN THE PAPER) -----------------------

TAB3 <- dat
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

# get the prevalence baseline -- cluster- and zone-specific  
TAB3_agg_baseline <- TAB3 %>% 
  filter(survey_round == 1) %>%
  group_by(cluster_number, distance_hotspot_categories) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE)) %>%
  rename(PCR_18S_result_BASELINE = PCR_18S_result) %>%
  as.data.frame()
TAB3 <- TAB3 %>% left_join(TAB3_agg_baseline, by = c("distance_hotspot_categories", "cluster_number"))
# make empty vector of predicted values
TAB3$PCR_18S_result_PRED <- NA

# do logistic regression to predict the outcome -- response time-specific 
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

TAB3_agg <- TAB3 %>% 
  group_by(cluster_number, distance_hotspot_categories, cluster_arm, survey_round) %>%
  summarise(PCR_18S_result = mean(PCR_18S_result, na.rm = TRUE),
            PCR_18S_result_PRED = mean(PCR_18S_result_PRED, na.rm = TRUE)) %>%
  mutate(PCR_18S_result_RESID = PCR_18S_result - PCR_18S_result_PRED) %>%
  as.data.frame()
dim(TAB3_agg)



# subset data from to keep survey_round == 2 only
## vals unadjusted 
TAB1_agg$y <- TAB1_agg$PCR_18S_result 
TAB1_sr2_dhc0 <- TAB1_agg %>% filter(survey_round == 2, distance_hotspot_categories == 0)  
TAB1_sr2_dhc1 <- TAB1_agg %>% filter(survey_round == 2, distance_hotspot_categories == 1)  
TAB1_sr2_dhc2 <- TAB1_agg %>% filter(survey_round == 2, distance_hotspot_categories == 2)  
## vals adjusted 
TAB3_agg$y <- TAB3_agg$PCR_18S_result_RESID 
TAB3_sr2_dhc0 <- TAB3_agg %>% filter(survey_round == 2, distance_hotspot_categories == 0)  
TAB3_sr2_dhc1 <- TAB3_agg %>% filter(survey_round == 2, distance_hotspot_categories == 1)  
TAB3_sr2_dhc2 <- TAB3_agg %>% filter(survey_round == 2, distance_hotspot_categories == 2)  
# list of data frame
dat_list <- list(
  TAB1_sr2_dhc0,
  TAB1_sr2_dhc1,
  TAB1_sr2_dhc2,
  TAB3_sr2_dhc0,
  TAB3_sr2_dhc1,
  TAB3_sr2_dhc2
)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# DO UPSTRAP
dat_nrow  <- nrow(TAB1_sr2_dhc0)
B_boot    <- 1000 * 100
n1_min    <- 5
n1_max    <- 30
n1_grid     <- n1_min:n1_max
n1_grid_l   <- length(n1_grid)
n1_grid_max <- max(n1_grid)

out_power   <- numeric()
out_dhc     <- numeric()
out_yval    <- numeric()
out_n1      <- numeric()

# data index corresponding to control/intervention arms
dat_rowidx_ctrl <- which(cluster_arm == "CONTROL") # [1] 1 2 3 5 9
dat_rowidx_itrv <- which(cluster_arm == "INTERVENTION") # [1]  4  6  7  8 10

set.seed(123)
for (i in 1:n1_grid_l){ # i <- 5; j <-1
  t1 <- Sys.time()
  n1_i <- n1_grid[i]
  # sample row index for current upstrap sampling: control, intervention 
  boot_rowidx_ctrl <- matrix(sample(x = dat_rowidx_ctrl, size = (B_boot * n1_i), replace = TRUE), nrow = B_boot)
  boot_rowidx_itrv <- matrix(sample(x = dat_rowidx_itrv, size = (B_boot * n1_i), replace = TRUE), nrow = B_boot)
  # repeat procedure B_boot-times 
  out_B_boot <- matrix(NA, nrow = B_boot, ncol = 6)
  for (j in 1:B_boot){
    if (j %% 1000 == 0) print(j)
    ij1_idx <- boot_rowidx_ctrl[j, ]
    ij2_idx <- boot_rowidx_itrv[j, ]
    ij_idx  <- c(ij1_idx, ij2_idx) 
    out_B_boot[j, ] <- sapply(dat_list, function(dat_tmp){
      (t.test(y ~ cluster_arm, data = dat_tmp[ij_idx, ]))$p.value < 0.05
    })
  }
  # compute estimated power 
  out_power_i <- colMeans(out_B_boot)
  # store simulation results
  out_power <- c(out_power, out_power_i)
  out_dhc  <- c(out_dhc, c(0:2, 0:2))
  out_yval <- c(out_yval, c(rep("x_weeks", 3), rep("x_weeks_adj", 3)))
  out_n1   <- c(out_n1, rep(n1_i, 6))
  t2 <- Sys.time()
  t_diff <- as.numeric(t2 - t1, unit = "secs")
  message(paste0("i: ", i, " | Time elapsed: ", round(t_diff), " s "))
}


out_power_df <- data.frame(
  out_n1 = out_n1, 
  out_power = out_power, 
  out_yval = out_yval,
  out_dhc = out_dhc)
out_power_df

out_power_df_fpath <- paste0(here::here(), "/use_case_examples/results/2020-12-14-upstrap_bousema2016.rds")
saveRDS(object = out_power_df, file = out_power_df_fpath)








