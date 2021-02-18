
#' This script is used to determine, roughly, grids of N0 sample size such that 
#' it is 10, 30, 50, 70, 90% power in the following problems: 
#' - one-sample t-test
#' - two-sample t-test
#' 

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(lme4)
library(lmerTest)
library(geepack)

# one-sample t-test ------------------------------------------------------------

x = power_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95)
y = sapply(power_grid, function(power_tmp) 
  power.t.test(power = power_tmp, delta = 0.3, sd = 1, type = "one.sample")$n)
data.frame(x, y = round(y))
# x   y
# 1 0.10   7
# 2 0.30  25
# 3 0.50  45
# 4 0.70  71
# 5 0.90 119
# 6 0.95 146



# two-sample t-test ------------------------------------------------------------

x = power_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95)
y = sapply(power_grid, function(power_tmp) 
  power.t.test(power = power_tmp, delta = 0.3, sd = 1, type = "two.sample")$n)
data.frame(x, y = round(y))
# x   y
# 1 0.10  11
# 2 0.30  47
# 3 0.50  86
# 4 0.70 138
# 5 0.90 234
# 6 0.95 290



# LMM (1 + 0 cov) --------------------------------------------------------------

coef_x1 <- 0.5
tau2    <- 1
sigma2  <- 1
N       <- 50   # sample size of each of the two arms
ni      <- 3
N1_min  <- 10
N1_max  <- 150

# define N1 grid
N1_grid   <- seq(from = N1_min, to = N1_max, by = 1)
N1_grid_l <- length(N1_grid)

rep_R <- 1000

# TMP
mat_out_LMM  <- matrix(NA, nrow = rep_R, ncol = N1_grid_l)
mat_out_GEE  <- matrix(NA, nrow = rep_R, ncol = N1_grid_l)

set.seed(123)
t1 <- Sys.time()
for (rep_i in 1:rep_R){
  print(paste0("r_rep: ", rep_i))

  # deterministic quantities
  N_tmp <- N1_max
  subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
  subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
  x1_i      <- c(rep(1, N_tmp), rep(0, N_tmp))
  # simulated quantities
  b0_i      <- rnorm(n = (N_tmp * 2), mean = 0, sd = tau2)
  eps_ij    <- rnorm(n = (N_tmp * 2 * ni), mean = 0, sd = sigma2)
  # simulated/generated data frame variables 
  subjid_ij     <- rep(subjid_i, each = ni) 
  subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
  x1_ij     <- rep(x1_i, each = ni) 
  b0_ij     <- rep(b0_i, each = ni) 
  y_ij      <- b0_ij + coef_x1 * x1_ij + eps_ij
  dat       <- data.frame(y = y_ij, x1 = x1_ij, subjid = subjid_ij, subjid_arm = subjid_arm_ij)
  
  # iterate over number of independent units (subjects) in the same arm
  for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 10
    N1_tmp <- N1_grid[N1_grid_idx]
    dat_sub <- dat %>% filter(subjid_arm <= N1_tmp)
    # # LMM
    # fit_LMM   <- lmer(y ~ x1 + (1 | subjid), data = dat_sub)
    # fit_LMM_s <- summary(fit_LMM)
    # fit_LMM_coefpval  <- fit_LMM_s$coefficients[2, 5]
    # mat_out_LMM[rep_i, N1_grid_idx] <- (fit_LMM_coefpval < 0.05) * 1
    # GEE
    fit_GEE_formula <- formula(y ~ x1)
    fit_GEE <- geeglm(formula = fit_GEE_formula, family = gaussian(link = "identity"), 
                      data = dat_sub, id = dat_sub$subjid, corstr = "exchangeable")
    fit_GEE_s <- summary(fit_GEE)
    fit_GEE_coefest_pval  <- fit_GEE_s$coefficients[2, 4]
    mat_out_GEE[rep_i, N1_grid_idx] <- (fit_GEE_coefest_pval < 0.05) * 1
  }
}
t2 <- Sys.time()
t2 - t1
# Time difference of 45.60784 mins

est_power_df <- data.frame(
  N1 = N1_grid,
  est_power_GEE = matrixStats::colMeans2(mat_out_GEE)
)
est_power_df

# target df 
target_power_df <- data.frame(
  target_power = c(0.3, 0.5, 0.7, 0.9, 0.95)
)
est_power_df %>% 
  full_join(target_power_df, by = character()) %>%
  mutate(power_diff = abs(est_power_GEE - target_power)) %>%
  group_by(target_power) %>%
  filter(power_diff == min(power_diff)) %>%
  summarise(N1 = mean(N1))
# target_power    N1
# <dbl> <dbl>
# 1         0.3    20 
# 2         0.5    41 
# 3         0.7    62 
# 4         0.9   111 
# 5         0.95  130.



# GLMM (1 + 0 cov) --------------------------------------------------------------

# define experiment params 
coef_x1 <- 0.5
tau2    <- 1
sigma2  <- 1
N       <- 50   # sample size of each of the two arms
ni      <- 3
N1_min  <- 50
N1_max  <- 150

# define N1 grid
N1_grid   <- seq(from = N1_min, to = N1_max, by = 3)
N1_grid_l <- length(N1_grid)

# deterministic quantities
N_tmp <- N1_max
subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
x1_i         <- c(rep(1, N_tmp), rep(0, N_tmp))


# step 1: choose the intercept ----------------------------------------------------

R_rep <- 1000
coef_0_grid <- seq(-1, 1, length.out = 100)
out_mat <- matrix(NA, nrow = R_rep, ncol = length(coef_0_grid))

for (rep_i in 1 : R_rep){
  for (coef_idx in 1 : length(coef_0_grid)){
    coef_0 <- coef_0_grid[coef_idx]
    # select constant to get ~ 50% 1's
    b0_i      <- rnorm(n = (N_tmp * 2), mean = 0, sd = tau2)
    # eps_ij    <- rnorm(n = (N_tmp * 2 * ni), mean = 0, sd = sigma2)
    # simulated/generated data frame variables 
    subjid_ij     <- rep(subjid_i, each = ni) 
    subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
    x1_ij     <- rep(x1_i, each = ni) 
    b0_ij     <- rep(b0_i, each = ni) 
    XB_ij     <- coef_0 + b0_ij + coef_x1 * x1_ij
    p_ij      <- 1/(1 + exp(-XB_ij))
    y_ij      <- rbinom(n = length(p_ij), size = 1, prob = p_ij)
    out_mat[rep_i, coef_idx] <- mean(y_ij)
  }
}
out_vec <- matrixStats::colMeans2(out_mat)
coef_0_grid[which.min(abs(out_vec - 0.5))]
# [1] -0.2525253


# step 2: N0, N1 grid ----------------------------------------------------0-----

rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(lme4)
library(lmerTest)
library(geepack)

# define experiment params 
coef_0  <- -0.25
coef_x1 <- 0.5
tau2    <- 1
sigma2  <- 1
N       <- 50   # sample size of each of the two arms
ni      <- 3
N1_min  <- 10
N1_max  <- 250

# define N1 grid
N1_grid   <- seq(from = N1_min, to = N1_max, by = 1)
N1_grid_l <- length(N1_grid)

# deterministic quantities
N_tmp <- N1_max
subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
x1_i         <- c(rep(1, N_tmp), rep(0, N_tmp))


rep_R <- 1000

# mat_out_GLMM <- matrix(NA, nrow = rep_R, ncol = N1_grid_l)
mat_out_GEE  <- matrix(NA, nrow = rep_R, ncol = N1_grid_l)

t1 <- Sys.time()
set.seed(123)
for (rep_i in 1 : rep_R){
  message(paste0("rep_i: ", rep_i))
  
  # simulated quantities
  b0_i      <- rnorm(n = (N_tmp * 2), mean = 0, sd = tau2)
  # simulated/generated data frame variables 
  subjid_ij     <- rep(subjid_i, each = ni) 
  subjid_arm_ij <- rep(subjid_arm_i, each = ni) 
  x1_ij     <- rep(x1_i, each = ni) 
  b0_ij     <- rep(b0_i, each = ni) 
  XB_ij     <- coef_0 + b0_ij + coef_x1 * x1_ij
  p_ij      <- 1/(1 + exp(-XB_ij))
  y_ij      <- rbinom(n = length(p_ij), size = 1, prob = p_ij)
  dat       <- data.frame(y = y_ij, x1 = x1_ij, subjid = subjid_ij, subjid_arm = subjid_arm_ij)
  
  for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 10
    # message(paste0("rep_i: ", rep_i, ", N1_grid_idx: ", N1_grid_idx))
    
    N1_tmp <- N1_grid[N1_grid_idx]
    dat_sub <- dat %>% filter(subjid_arm <= N1_tmp)
    
    # # GLMM
    # fit_GLMM   <- glmer(y ~ x1 + (1 | subjid), data = dat_sub, family = binomial)
    # fit_GLMM_s <- summary(fit_GLMM)
    # fit_GLMM_coefest_pval  <- fit_GLMM_s$coefficients[2, 4]
    # mat_out_GLMM[rep_i, N1_grid_idx] <- (fit_GLMM_coefest_pval < 0.05) * 1
    
    # GEE
    fit_GEE_formula <- formula(y ~ x1)
    fit_GEE <- geeglm(formula = fit_GEE_formula, family = binomial(link = "logit"), 
                      data = dat_sub, id = dat_sub$subjid, corstr = "exchangeable")
    fit_GEE_s <- summary(fit_GEE)
    fit_GEE_coefest_pval  <- fit_GEE_s$coefficients[2, 4]
    mat_out_GEE[rep_i, N1_grid_idx] <- (fit_GEE_coefest_pval < 0.05) * 1
  }
}
t2 <- Sys.time()
t2 - t1

est_power_df <- data.frame(
  N1 = N1_grid,
  est_power_GEE = matrixStats::colMeans2(mat_out_GEE)
)
est_power_df

# target df 
target_power_df <- data.frame(
  target_power = c(0.3, 0.5, 0.7, 0.9)
)
est_power_df %>% 
  full_join(target_power_df, by = character()) %>%
  mutate(power_diff = abs(est_power_GEE - target_power)) %>%
  group_by(target_power) %>%
  filter(power_diff == min(power_diff)) %>%
  summarise(N1 = mean(N1))

