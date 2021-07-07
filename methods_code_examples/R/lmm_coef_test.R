rm(list = ls())

# simulation parameters
N       <- 40   
ni      <- 3
coef_x1 <- 1
coef_x2 <- 0.01
tau2    <- 1
sigma2  <- 1
R_boot  <- 1000
# simulate sample
set.seed(1)
subjid_i      <- 1:N                       # subject ID unique in data set 
x1_i          <- rep(c(0, 1), times = N/2) # indicator of a "trial arm"
x2_i          <- runif(n = N, min = 18, max = 100)
b0_i          <- rnorm(n = N, mean = 0, sd = tau2)
eps_ij        <- rnorm(n = N * ni, mean = 0, sd = sigma2)
subjid_ij     <- rep(subjid_i, each = ni) 
x1_ij         <- rep(x1_i, each = ni) 
x2_ij         <- rep(x2_i, each = ni) 
b0_ij         <- rep(b0_i, each = ni) 
y_ij          <- b0_ij + coef_x1 * x1_ij + eps_ij
dat           <- data.frame(y = y_ij, x1 = x1_ij, x2 = x2_ij, subjid = subjid_ij)
head(dat, 4)
#           y x1       x2 subjid
# 1 3.3205951  0 39.77171      1
# 2 0.8797374  0 39.77171      1
# 3 1.6087167  0 39.77171      1
# 4 1.8101385  1 48.51416      2

# Get observed effect size (x1 covariate coefficient estimate)
library(lme4)
library(lmerTest)
fit  <- lmer(y ~ x1 + x2 + (1 | subjid), data = dat)
fixef(fit)["x1"]
# 0.6890758 


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
  fit_rr  <- lmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr)
  pval_rr <- summary(fit_rr)$coef["x1", 5]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.536


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: 60, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat$subjid), size = 60, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat[dat$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : 60, each = ni)
  fit_rr  <- lmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr)
  pval_rr <- summary(fit_rr)$coef["x1", 5]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.714


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# target sample size: 60, target effect size: 1.2
dat_upd <- dat
dat_upd$y <- dat_upd$y + (1.2 - fixef(fit)["x1"]) * dat_upd$x1
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  print(rr)
  dat_rr_subj_id <- sample(unique(dat_upd$subjid), size = 60, replace = TRUE)
  dat_rr  <- lapply(dat_rr_subj_id, function(subjid_tmp) dat_upd[dat_upd$subjid == subjid_tmp, ]) 
  dat_rr  <- do.call("rbind", dat_rr)
  # make new subject ID so as to treat subjects resampled >1 as unique ones
  dat_rr$subjid <- rep(1 : 60, each = ni)
  fit_rr  <- lmer(y ~ x1 + x2 + (1 | subjid), data = dat_rr)
  pval_rr <- summary(fit_rr)$coef["x1", 5]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.994

