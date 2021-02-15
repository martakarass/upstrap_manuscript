

rm(list = ls())
options(scipen = 999)
library(data.table)

mu         <- 0.1 
n0         <- 200
rep_n      <- 1000
rep_idx_1_grid  <- 1 + seq(0, by = rep_n, length.out = 100)
B_boot     <- 10000
n_grid_min <- n0
n_grid_max <- 1200

vec <- numeric()
for (rep_idx_1 in rep_idx_1_grid){
  vec_i <- paste0(
    "MU_", mu, 
    "_N0_", n0, 
    "_R_", rep_n, 
    "_BBOOT_", B_boot,
    "_N1MIN_", n_grid_min,
    "_N1MAX_", n_grid_max,
    "_REPIDX1_", rep_idx_1)
  vec <- c(vec, vec_i)
}

df <- data.frame(vec = vec)
df_fname <- "2020-11-23-ttest_normal_mean_params.txt"
df_fpath <- paste0("/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/numerical_experiments/utils_CL/", df_fname)
df

fwrite(as.data.table(df), df_fpath, col.names = FALSE)
