rm(list = ls())

# simulation parameters
N <- 30
R_boot <- 1000 

# simulate sample
set.seed(1)
x <- rnorm(n = N, mean = 0.3, sd = 1)
mean(x)
# [1] 0.3824582


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: as observed, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.593

power.t.test(n = N, delta = mean(x), sd = sd(x), type = "one.sample")$power
# [1] 0.5914924


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: M = 40, target effect size: as observed
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x, size = 40, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.716

power.t.test(n = 40, delta = mean(x), sd = sd(x), type = "one.sample")$power
# [1] 0.7232818


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# estimate power 
# target sample size: M = 40, target effect size: 0.5
x_upd <- x + (0.5 - mean(x)) * 1
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x_upd, size = 40, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.917

power.t.test(n = 40, delta = 0.5, sd = sd(x), type = "one.sample")$power
# [1] 0.9156754


