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

# Get observed effect size (x1 covariate coefficient estimate)
library(lme4)
library(lmerTest)
fit  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat, family = binomial)
fixef(fit)["x1"]
# 0.8574488


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: as observed, target effect size: as observed
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
# target sample size: 100, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat$subjid), size = 100, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat[dat$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : 100, each = ni)
  fit_rr  <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)
  pval_rr <- summary(fit_rr)$coef["x1", 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.794


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: 100, target effect size: 1.2

dat_upd <- dat
# set re.form = NULL to condition on all random effects when "predicting" the link response
dat_upd$link_orig <- predict(fit, type = "link", re.form = NULL) 
# update link function to express the assumed target effect size 
dat_upd$link_upd  <- dat_upd$link_orig + (1.2 - fixef(fit)["x1"]) * dat_upd$x1
# update fitted response (probability of Y=1)
dat_upd$res_upd   <- 1/(1 + exp(-dat_upd$link_upd))
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat_upd$subjid), size = 100, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat_upd[dat_upd$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : 100, each = ni)
  # simulate response on resampled data, assuming target effect size  
  dat_rr$y <- rbinom(n = nrow(dat_rr), size = 1, prob = dat_rr$res_upd)
  fit_rr   <- glmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr, family = binomial)
  pval_rr  <- summary(fit_rr)$coef["x1", 4]
  out[rr]  <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.992


