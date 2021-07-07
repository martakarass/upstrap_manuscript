rm(list = ls())

# simulation parameters
N       <- 80
coef_x0 <- -0.2
coef_x1 <- 0.5
coef_x2 <- 0.1 
coef_x3 <- -0.01 
sigma2  <- 1
R_boot  <- 1000 
# simulate sample
set.seed(1)
subjid_i <- 1 : N    
x1_i     <- rep(c(0, 1), times = N/2)
x2_i     <- rbinom(n = N, size = 1, prob = 0.5)
x3_i     <- runif(n = N, min = 18, max = 100)
XB_i     <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) + (coef_x3 * x3_i)
p_i      <- 1/(1 + exp(-XB_i))
y_i      <- rbinom(n = length(p_i), size = 1, prob = p_i)
dat      <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i, subjid = subjid_i)
head(dat, 3)
#   y x1 x2       x3 subjid
# 1 0  0  0 53.64208      1
# 2 0  1  0 76.42620      2
# 3 0  0  1 50.79954      3

# Get observed effect size (x1 covariate coefficient estimate)
fit  <- glm(y ~ x1 + x2 + x3, data = dat, family = binomial(link = "logit"))
coef(fit)["x1"]
# 0.7354476 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: as observed, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat), replace = TRUE)
  dat_rr  <- dat[dat_rr_idx, ]
  fit_rr  <- glm(y ~ x1 + x2 + x3, data = dat_rr, 
                 family = binomial(link = "logit"))
  pval_rr <- summary(fit_rr)$coef[2, 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.356


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: 120, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat), size = 120, replace = TRUE)
  dat_rr  <- dat[dat_rr_idx, ]
  fit_rr  <- glm(y ~ x1 + x2 + x3, data = dat_rr, family = binomial(link = "logit"))
  pval_rr <- summary(fit_rr)$coef[2, 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.441


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: 80, target effect size: 1.0
dat_upd <- dat
dat_upd$link_orig <- predict(fit, type = "link") 
dat_upd$link_upd  <- dat_upd$link_orig  + (1 - coef(fit)["x1"]) * dat_upd$x1
dat_upd$res_upd   <- 1/(1 + exp(-dat_upd$link_upd))
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat_upd), size = 80, replace = TRUE)
  dat_rr   <- dat_upd[dat_rr_idx, ]
  # simulate response on resampled data, assuming target effect size  
  dat_rr$y <- rbinom(n = nrow(dat_rr), size = 1, prob = dat_rr$res_upd)
  fit_rr   <- glm(y ~ x1 + x2 + x3, data = dat_rr, family = binomial(link = "logit"))
  pval_rr  <- summary(fit_rr)$coef[2, 4]
  out[rr]  <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.565
