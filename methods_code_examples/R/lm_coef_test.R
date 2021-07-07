rm(list = ls())

# simulation parameters
N       <- 40
coef_x0 <- 0 
coef_x1 <- 0.6
coef_x2 <- 0.3 
coef_x3 <- -0.1
sigma2  <- 1
R_boot  <- 1000 
# simulate sample
set.seed(123)
subjid_i  <- 1 : N   # subject ID unique in data set 
x1_i      <- rbinom(n = N, size = 1, prob = 0.5)
x2_i      <- rbinom(n = N, size = 1, prob = 0.5)
x3_i      <- runif(n = N, min = 18, max = 100)
eps_i     <- rnorm(N, sd = sqrt(sigma2))
y_i       <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) + 
  (coef_x3 * x3_i) + eps_i
dat       <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i,
                        subjid = subjid_i)
head(dat, 3)
# y           x1 x2       x3 subjid
# 1 -2.662590  0  1 53.64208      1
# 2 -7.381860  0  1 76.42620      2
# 3 -3.590214  1  1 50.79954      3

# Get observed effect size (x1 covariate coefficient estimate)
fit  <- lm(y ~ x1 + x2 + x3, data = dat)
coef(fit)["x1"]
# 0.3949275 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: as observed, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat), replace = TRUE)
  dat_rr  <- dat[dat_rr_idx, ]
  fit_rr  <- lm(y ~ x1 + x2 + x3, data = dat_rr)
  pval_rr <- summary(fit_rr)$coef[2, 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.286


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# target sample size: 80, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat), size = 80, replace = TRUE)
  dat_rr  <- dat[dat_rr_idx, ]
  fit_rr  <- lm(y ~ x1 + x2 + x3, data = dat_rr)
  pval_rr <- summary(fit_rr)$coef["x1", 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.511


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# target sample size: 80, target effect size: 0.5 
dat_upd <- dat
dat_upd$y <- dat_upd$y + (0.5 - coef(fit)["x1"]) * dat_upd$x1
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  dat_rr_idx <- sample(1 : nrow(dat_upd), size = 80, replace = TRUE)
  dat_rr  <- dat_upd[dat_rr_idx, ]
  fit_rr  <- lm(y ~ x1 + x2 + x3, data = dat_rr)
  pval_rr <- summary(fit_rr)$coef[2, 4]
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.686

