rm(list = ls())

# simulation parameters
N <- 30 
R_boot <- 1000 
# simulate sample
set.seed(1)
x1 <- rnorm(N, mean = 0, sd = 1)
x2 <- rnorm(N, mean = 0.3, sd = 1) 

mean(x2) - mean(x1)
# [1] 0.3503164


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: as observed, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x1_rr <- sample(x1, replace = TRUE)
  x2_rr <- sample(x2, replace = TRUE)
  pval_rr <- t.test(x1_rr, x2_rr, paired = FALSE, var.equal = TRUE)$p.value
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.349

obs_delta <- mean(x2) - mean(x1)
obs_var_pooled <- (var(x1) * (N - 1) + var(x2) * (N - 1)) / (N + N - 2)
obs_sd_pooled <- sqrt(obs_var_pooled)
power.t.test(n = N, delta = obs_delta, sd = obs_sd_pooled, type = "two.sample")$power
# [1] 0.3400754


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: M = 40, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x1_rr <- sample(x1, size = 40, replace = TRUE)
  x2_rr <- sample(x2, size = 40, replace = TRUE)
  pval_rr <- t.test(x1_rr, x2_rr, paired = FALSE, var.equal = TRUE)$p.value
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.44

power.t.test(n = 40, delta = obs_delta, sd = obs_sd_pooled, type = "two.sample")$power
# [1] 0.4344143


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: M = 40, target effect size: 0.5
x1_upd <- x1
x2_upd <- x2 + (0.5 - (mean(x2) - mean(x1))) * 1
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x1_rr <- sample(x1_upd, size = 40, replace = TRUE)
  x2_rr <- sample(x2_upd, size = 40, replace = TRUE)
  pval_rr <- t.test(x1_rr, x2_rr, paired = FALSE, var.equal = TRUE)$p.value
  out[rr] <- (pval_rr < 0.05)
}
mean(out)
# [1] 0.752

power.t.test(n = 40, delta = 0.5, sd = obs_sd_pooled, type = "two.sample")$power
# [1] 0.7262939


