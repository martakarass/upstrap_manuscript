rm(list = ls())

# simulation parameters
N       <- 80   
ni      <- 3
coef_x1 <- 0.8
coef_x2 <- 0.01
tau2    <- 1
sigma2  <- 1
R_boot  <- 1000
# simulate sample
set.seed(1)
subjid_i  <- 1:N                       # subject ID unique in data set 
x1_i      <- rep(c(0, 1), times = N/2) # indicator of a "trial arm"
x2_i      <- runif(n = N, min = 18, max = 100)
b0_i      <- rnorm(n = N, mean = 0, sd = tau2)
subjid_ij <- rep(subjid_i, each = ni) 
x1_ij     <- rep(x1_i, each = ni) 
x2_ij     <- rep(x2_i, each = ni) 
b0_ij     <- rep(b0_i, each = ni) 
XB_ij     <- b0_ij + coef_x1 * x1_ij 
p_ij      <- 1/(1 + exp(-XB_ij))
y_ij      <- rbinom(n = length(p_ij), size = 1, prob = p_ij)
dat       <- data.frame(y = y_ij, x1 = x1_ij, subjid = subjid_ij)
dat       <- data.frame(y = y_ij, x1 = x1_ij, x2 = x2_ij, subjid = subjid_ij)
head(dat, 4)
#   y x1       x2 subjid
# 1 0  0 39.77171      1
# 2 1  0 39.77171      1
# 3 1  0 39.77171      1
# 4 1  1 48.51416      2


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(lme4)
library(lmerTest)
fit  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat, family = binomial)
fixef(fit)["x1"]
# 0.8574488


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat$subjid), replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat[dat$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : N, each = ni)
  fit_rr  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)
  pval_rr <- summary(fit_rr)$coef["x1", 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# 0.652


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

M <- 100
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat$subjid), size = M, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat[dat$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : M, each = ni)
  fit_rr  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)
  pval_rr <- summary(fit_rr)$coef["x1", 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.794


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

M <- 100
set.seed(1)
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat$subjid), size = M, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat[dat$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : M, each = ni)
  # fit model, update coefficient for the size effect of interest 
  fit_rr_sim_y  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)
  fixef(fit_rr_sim_y)["x1"] <- 1.2
  # simulate new outcome 
  dat_rr$y <- simulate(fit_rr_sim_y, nsim = 1, newdata = dat_rr)[[1]]
  # fit "final" model 
  fit_rr  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)  
  pval_rr <- summary(fit_rr)$coef["x1", 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.957
