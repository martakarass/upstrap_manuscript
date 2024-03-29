rm(list = ls())

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


N <- 30
R_boot <- 1000 
set.seed(1)
x <- rnorm(n = N, mean = 0.3, sd = 1)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

mean(x)
# [1] 0.3824582

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.593


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


M <- 40
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x, size = M, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.716


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# update sample x to represent the target effect size
M <- 40
x_upd <- x + (0.5 - mean(x)) * 1
out <- rep(NA, R_boot)
for (rr in 1 : R_boot){
  x_rr <- sample(x_upd, size = M, replace = TRUE)
  out[rr] <- (t.test(x_rr)$p.value < 0.05)
}
mean(out)
# [1] 0.917



