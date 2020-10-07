

rm(list = ls())

options(scipen = 999)

library(data.table)

R      <- 100000
BBOOT   <- 10000
N0_grid <- c(100, 200, 500, 1000, 2000)
N1_frac_grid <- c(0.9, 0.8, 0.4, 0.2, 0.1)

vec <- numeric()
for (N0 in N0_grid){
  for(N1_frac in N1_frac_grid){
    N1 <- round(N0 * N1_frac)
    vec_i <- paste0("N0_", N0, "_N1_", N1, "_R_", R, "_BBOOT_", BBOOT)
    vec <- c(vec, vec_i)
  }
}

df <- data.frame(vec = vec)
df_fpath <- "/Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/numerical_experiments/utils_CL/2020-10-06-N0_N1_R_BOOT_params.txt"
df

fwrite(as.data.table(df), df_fpath, col.names = FALSE)
