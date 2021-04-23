
rm(list = ls())
library(data.table)
library(dplyr)

dat_dir <- paste0(here::here(), "/use_case_examples/data_processed/bousema2016")

# list all files in data dir 
dat_fnames <- list.files(dat_dir, full.names = TRUE)

dat_baseline   <- fread(dat_fnames[1]) %>% as.data.frame()
dat_interv_cov <- fread(dat_fnames[2]) %>% as.data.frame()
dat_interv     <- fread(dat_fnames[3]) %>% as.data.frame()
dat_entomology <- fread(dat_fnames[4]) %>% as.data.frame()
names_list <- list(
  names(dat_baseline),
  names(dat_interv_cov),
  names(dat_interv),
  names(dat_entomology)
)
sapply(names_list, function(df_names){ "PCR_18S_result" %in% df_names})
sapply(names_list, function(df_names){ "cluster_number" %in% df_names})
sapply(names_list, function(df_names){ "hotspot_evaluation" %in% df_names})
sapply(names_list, function(df_names){ "distance_hotspot_categories" %in% df_names})

# data which has all column names used 
TAB1 <- dat_interv

# RESULT OF PCR / PCR_18S_result
summary(as.factor(TAB1$PCR_18S_result))
# PER CLUSTER : cluster_number
summary(as.factor(TAB1$cluster_number))
# PER HOTSPOT/EVALUATION : hotspot_evaluation
summary(as.factor(TAB1$hotspot_evaluation))
# PER DISTANCE FROM HOTSPOT : distance_hotspot_categories
summary(as.factor(TAB1$distance_hotspot_categories))



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----       left hand side of the Naudet table cell   -----------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

# PREVALENCE CUMULATIVE (?) (8 WEEKS and 16 WEEKS) -----------------------------

# @MK: following the analysis done in the paper 
# PREVALENCE IN EACH HOTSPOT (8 WEEKS and 16 WEEKS)
CLUSTER <- 0
PREVALENCE8 <- 0
PREVALENCE16 <- 0

TABPREVALENCE <- data.frame(CLUSTER, PREVALENCE8, PREVALENCE16)
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot" & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot"& TAB1$survey_round==2]))
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot" & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot"& TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TABPREVALENCE <- TABPREVALENCE[-1,]


# INTERVENTION; CONTROL
# add intervention/control  label to the cluster label
INTERVENTION <- c("CONTROL","CONTROL","CONTROL",
                  "INTERVENTION","CONTROL","INTERVENTION",
                  "INTERVENTION","INTERVENTION","CONTROL",
                  "INTERVENTION")
TABPREVALENCE <- cbind(TABPREVALENCE,INTERVENTION)


# SUMMARY AS SUGGESTED IN THE PAPER # t.test 8 WEEKS : 
# CLOSELY REPRODUCED  pVALUES AND 95 CONFINDENCE INTERVALS ARE SLIGTHLY DIFFERENTS
t.test(TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
# -0.02536956  0.22948491 ; p-value = 0.09488

t.test(TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
#  -0.05284218  0.13872766; p-value = 0.328


## -----------------------------------------------------------------------------
# PREVALENCE IN EACH EVALUATION ZONE 1 (8 WEEKS and 16 WEEKS) ------------------

#' dichotomized into two categories based on distance from the hotspot boundary (1–249 m and 250–500 m).
#' 
#' hotspot
#' evaluation zone 1 (TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==1)
#' evaluation zone 2 (TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==2)


CLUSTER <- 0
PREVALENCE8 <- 0
PREVALENCE16 <- 0
TABPREVALENCE <- data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==1 & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==1 & TAB1$survey_round==2])) 
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==1 & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==1 & TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TABPREVALENCE <- TABPREVALENCE[-1,]
TABPREVALENCE <- cbind(TABPREVALENCE, INTERVENTION)


t.test(TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
# -0.03302027  0.10439845; p-value = 0.2394

t.test(TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
# -0.07593394  0.09657864, p-value = 0.7778


## -----------------------------------------------------------------------------
# PREVALENCE IN EACH EVALUATION ZONE 2 (8 WEEKS and 16 WEEKS) ------------------

CLUSTER <- 0
PREVALENCE8 <- 0
PREVALENCE16 <- 0

TABPREVALENCE <- data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==2 & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==2 & TAB1$survey_round==2])) 
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==2 & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==2 & TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TABPREVALENCE <- TABPREVALENCE[-1,]
TABPREVALENCE <- cbind(TABPREVALENCE,INTERVENTION)


t.test(TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE8[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
# -0.02132842  0.10165905, p-value = 0.1698

t.test(TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="CONTROL"],
       TABPREVALENCE$PREVALENCE16[TABPREVALENCE$INTERVENTION=="INTERVENTION"])
# -0.08387298  0.10468770, p-value = 0.8046



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# reproduce the t-test results

INTERVENTION <- c("CONTROL","CONTROL","CONTROL",
                  "INTERVENTION","CONTROL","INTERVENTION",
                  "INTERVENTION","INTERVENTION","CONTROL",
                  "INTERVENTION")

# TAB1: hotspot
TABPREVALENCE <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot" & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot"& TAB1$survey_round==2]))
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot" & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Hotspot"& TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TAB_A <- cbind(TABPREVALENCE,INTERVENTION) 
TAB_A_ITRV <- TAB_A[TAB_A$INTERVENTION == "INTERVENTION", ]
TAB_A_CTRL <- TAB_A[TAB_A$INTERVENTION == "CONTROL", ]

# TAB2: eval_zone_1
TABPREVALENCE <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==1 & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==1 & TAB1$survey_round==2])) 
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==1 & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==1 & TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TAB_B <- cbind(TABPREVALENCE, INTERVENTION)
TAB_B_ITRV <- TAB_B[TAB_B$INTERVENTION == "INTERVENTION", ]
TAB_B_CTRL <- TAB_B[TAB_B$INTERVENTION == "CONTROL", ]

# TAB3: eval_zone_2
TABPREVALENCE <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (i in 1:10){
  PREVALENCE8 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==2 & TAB1$survey_round==2]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==2 & TAB1$survey_round==2])) 
  PREVALENCE16 <- sum(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation" & TAB1$distance_hotspot_categories==2 & TAB1$survey_round==3]))/length(na.omit(TAB1$PCR_18S_result[TAB1$cluster_number==i & TAB1$hotspot_evaluation=="Evaluation"& TAB1$distance_hotspot_categories==2 & TAB1$survey_round==3]))
  CLUSTER <- i
  TABPREVALENCE <- rbind(TABPREVALENCE,data.frame(CLUSTER,PREVALENCE8,PREVALENCE16)) 
}
TAB_C <- cbind(TABPREVALENCE, INTERVENTION)
TAB_C_ITRV <- TAB_C[TAB_C$INTERVENTION == "INTERVENTION", ]
TAB_C_CTRL <- TAB_C[TAB_C$INTERVENTION == "CONTROL", ]


B_boot    <- 10000
n1_min    <- 5
n1_max    <- 30
n1_grid     <- n1_min:n1_max
n1_grid_l   <- length(n1_grid)
n1_grid_max <- max(n1_grid)

# matrix to store upstrap results for each n1 value
out_power <- numeric()
out_zone  <- numeric()
out_n1_grid <- numeric()

set.seed(123)
for (i in 1:n1_grid_l){ # i <- 5; j <-1
  t1 <- Sys.time()
  # message(paste0("i: ", i))
  n1_i <- n1_grid[i]
  # sample row index: CONTROL, INTERVENTION
  boot_rowidx_i1 <- matrix(sample(x = 1:nrow(TAB_A_CTRL), size = (B_boot * n1_i), replace = TRUE), 
                          nrow = B_boot, ncol = n1_i)
  boot_rowidx_i2 <- matrix(sample(x = 1:nrow(TAB_A_ITRV), size = (B_boot * n1_i), replace = TRUE), 
                          nrow = B_boot, ncol = n1_i)
  # objects to store boot results
  out_B_boot_TAB_A <- rep(NA, B_boot)
  out_B_boot_TAB_B <- rep(NA, B_boot)
  out_B_boot_TAB_C <- rep(NA, B_boot)
  for (j in 1:B_boot){
    ij1_idx <- boot_rowidx_i1[j, ]
    ij2_idx <- boot_rowidx_i2[j, ]
    out_B_boot_TAB_A[j] <- (t.test(TAB_A_CTRL[ij1_idx, "PREVALENCE8"], TAB_A_ITRV[ij2_idx, "PREVALENCE8"]))$p.value < 0.05
    out_B_boot_TAB_B[j] <- (t.test(TAB_B_CTRL[ij1_idx, "PREVALENCE8"], TAB_B_ITRV[ij2_idx, "PREVALENCE8"]))$p.value < 0.05
    out_B_boot_TAB_C[j] <- (t.test(TAB_C_CTRL[ij1_idx, "PREVALENCE8"], TAB_C_ITRV[ij2_idx, "PREVALENCE8"]))$p.value < 0.05
  }
  out_power <- c(out_power, c(mean(out_B_boot_TAB_A), mean(out_B_boot_TAB_B), mean(out_B_boot_TAB_C)))
  out_zone <- c(out_zone, c("hotspot", "eval_zone_1", "eval_zone_2"))
  out_n1_grid <- c(out_n1_grid, rep(n1_i, 3))
  
  t2 <- Sys.time()
  t_diff <- as.numeric(t2 - t1, unit = "secs")
  message(paste0("i: ", i, " | Time elapsed: ", round(t_diff), " s "))
}


out_power_df <- data.frame(n1_grid = out_n1_grid, out_power = out_power, out_zone = out_zone)
out_power_df_fpath <- paste0(here::here(), "/use_case_examples/results/bousema2016_analysis1_outpower.rds")
saveRDS(object = out_power_df, file = out_power_df_fpath)











