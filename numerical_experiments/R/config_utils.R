
#' This script contains set of short functions utilized over multiple 
#' scripts. 

# function to get approximation of correction factor c4
# https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
get_c4_n <- function(n){
  c4_n <- 1 - (1/(4 * n)) -  (7/(32 * n^2)) - (19/(128 * n^3))
  return(c4_n)
}

# function to compute sd of "s" estimator 
get_sd_s <- function(n, sigma = 1){
  c4_n <- get_c4_n(n)
  sd_s <- sigma * sqrt(1 - c4_n^2)
  return(sd_s)
}

# function to compute E of "s" estimator 
get_mean_s <- function(n, sigma = 1){
  c4_n <- get_c4_n(n)
  mean_s <- sigma * c4_n
  return(mean_s)
}

message("The script config_utils.R was read.")
