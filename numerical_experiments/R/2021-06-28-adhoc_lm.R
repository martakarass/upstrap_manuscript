
rm(list = ls())
library(here)
library(tidyverse)
library(matrixStats)
library(simr)

# dir to save results 
res_fdir_agg  <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2021-06-28-lm_agg")
res_fdir_raw  <- paste0(here::here(), "/numerical_experiments/results_CL/2021-06-28-lm_raw")
res_fdir_plt  <- paste0(here::here(), "/numerical_experiments/results_figures/2021-06-28-lm")
# remove dirs if exist to make a room for new ones
unlink(res_fdir_agg, recursive = TRUE, force = TRUE)
unlink(res_fdir_raw, recursive = TRUE, force = TRUE)
unlink(res_fdir_plt, recursive = TRUE, force = TRUE)
# create dirs if any does not exist
dir.create(path = res_fdir_agg)
dir.create(path = res_fdir_raw)
dir.create(path = res_fdir_plt)
message(paste0("dir.exists(path = res_fdir_agg): ", dir.exists(path = res_fdir_agg)))
message(paste0("dir.exists(path = res_fdir_raw): ", dir.exists(path = res_fdir_raw)))
message(paste0("dir.exists(path = res_fdir_plt): ", dir.exists(path = res_fdir_plt)))

# experiment parameters
N0 <- 31
N1_min  <- 10
N1_max  <- 150
N1_grid   <- seq(from = N1_min, to = N1_max, by = 1)
N1_grid_l <- length(N1_grid)

coef_x0 <- 0 # set to zero 
coef_x1 <- 0.5
coef_x2 <- 0.3 # set to zero 
coef_x3 <- -0.1 # set to zero 
sigma2  <- 1

# number of repetitions of experiment 
R_rep   <- 200
# number of boot repetitions within one experiment, one setup
B_boot  <- 100

# effect of interest
effect1_grid <- c(0.3, 0.5, 0.7)

t1 <- Sys.time()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RUN SIMULATION

for (arrayjob_idx in 1 : R_rep){ # arrayjob_idx <- 1

    # seed seed for reproducibility of the results 
    set.seed(arrayjob_idx)
    message(paste0("arrayjob_idx: ", arrayjob_idx))
    # deterministic quantities
    N_tmp        <- N1_max
    subjid_i     <- 1:(N_tmp * 2)         # subject ID unique in whole data set 
    subjid_arm_i <- c(1:N_tmp, 1:N_tmp)   # subject ID unique in a trt arm
    x1_i         <- c(rep(1, N_tmp), rep(0, N_tmp))
    # simulated quantities
    x2_i         <- rbinom(n = N_tmp * 2, size = 1, prob = 0.5)
    x3_i         <- runif(n = N_tmp * 2, min = 18, max = 100)
    eps_i        <- rnorm(N_tmp * 2, sd = sqrt(sigma2))
    # simulated/generated data frame variables 
    y_i   <- coef_x0 + (coef_x1 * x1_i) + (coef_x2 * x2_i) + (coef_x3 * x3_i) + eps_i
    dat   <- data.frame(y = y_i, x1 = x1_i, x2 = x2_i, x3 = x3_i,
                        subjid = subjid_i, subjid_arm = subjid_arm_i)
    
    # make object to store simulation results 
    mat_out_all <- data.frame()
    
    # data with "observed" sample of size N0 
    datN0 <-  dat %>% filter(subjid_arm <= N0)
    datN0_1_idx <- which(datN0$x1 == 1)
    datN0_0_idx <- which(datN0$x1 == 0)
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH upstrap
    print("--- ESTIMATE POWER WITH upstrap")
    
    t_i_1 <- Sys.time()
    # define object to store values across B resamplings
    mat_out_boot <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
    # iterate over bootstrap repetitions 
    for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 10
        if (B_boot_idx %% 100 == 0){
            t_passed <- round(as.numeric(Sys.time() - t1, unit = "mins"))
            print(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", t_passed, " mins"))
        }
        # upsample data for current boot repetition (upstrap up to N1 max)
        dat_b <- rbind(
            datN0[sample(datN0_1_idx, size = N1_max, replace = TRUE), ], 
            datN0[sample(datN0_0_idx, size = N1_max, replace = TRUE), ])
        ## make new subj ID so as to treat resampled subjects as new ones
        dat_b$subjid <- 1 : (2 * N1_max)
        dat_b$subjid_arm <- c(1 : N1_max, 1 : N1_max)
        # iterate over N1 values grid
        for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 1
            tryCatch({
                fit   <- lm(y ~ x1 + x2 + x3, data = dat_b[dat_b$subjid_arm <= N1_grid[N1_grid_idx], ])
                fit_s <- summary(fit)
                fit_coefpval  <- fit_s$coefficients[2, 4]
                mat_out_boot[N1_grid_idx, B_boot_idx] <- (fit_coefpval < 0.05) * 1
            }, error = function(e) {message(e)})
        }
    }
    # add results to mat_out
    value <- rowMeans(mat_out_boot, na.rm = TRUE); rm(mat_out_boot)
    t_i_2 <- Sys.time()
    t_i_12_diff <- round(as.numeric(difftime(t_i_2, t_i_1, units = "secs")), 1)
    
    # store results to master file
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$name          <- "upstrap_power"
    mat_out_tmp$effect0       <- coef_x1
    mat_out_tmp$effect1       <- NA
    mat_out_tmp$value         <- value
    mat_out_tmp$elapsed_s     <- t_i_12_diff
    mat_out_all <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value, t_i_12_diff, t_i_1, t_i_2)
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH upstrap, modify effect 
    print("--- ESTIMATE POWER WITH upstrap, modify effect ")
    
    for (effect1 in effect1_grid){
    
        # update sample
        fit   <- lm(y ~ x1 + x2 + x3, data = datN0)
        effect_obs <- coef(fit)[2]
        datN0$y <- datN0$y + ((effect1 * datN0$x1) - (effect_obs * datN0$x1))
        
        datN0_1_idx <- which(datN0$x1 == 1)
        datN0_0_idx <- which(datN0$x1 == 0)
        
        # define object to store values across B resamplings
        mat_out_boot <- matrix(NA, nrow = N1_grid_l, ncol = B_boot)
        # iterate over bootstrap repetitions 
        for (B_boot_idx in 1:B_boot){ # B_boot_idx <- 10
            if (B_boot_idx %% 100 == 0){
                t_passed <- round(as.numeric(Sys.time() - t1, unit = "mins"))
                print(paste0("B_boot_idx: ", B_boot_idx, " [", round(B_boot_idx / B_boot * 100, 2), "%], ", t_passed, " mins"))
            }
            # upsample data for current boot repetition (upstrap up to N1 max)
            dat_b <- rbind(
                datN0[sample(datN0_1_idx, size = N1_max, replace = TRUE), ], 
                datN0[sample(datN0_0_idx, size = N1_max, replace = TRUE), ])
            ## make new subj ID so as to treat resampled subjects as new ones
            dat_b$subjid <- 1 : (2 * N1_max)
            dat_b$subjid_arm <- c(1 : N1_max, 1 : N1_max)
            # iterate over N1 values grid
            for (N1_grid_idx in 1 : N1_grid_l){ # N1_grid_idx <- 1
                tryCatch({
                    fit   <- lm(y ~ x1 + x2 + x3, data = dat_b[dat_b$subjid_arm <= N1_grid[N1_grid_idx], ])
                    fit_s <- summary(fit)
                    fit_coefpval  <- fit_s$coefficients[2, 4]
                    mat_out_boot[N1_grid_idx, B_boot_idx] <- (fit_coefpval < 0.05) * 1
                }, error = function(e) {message(e)})
            }
        }
        # add results to mat_out
        value <- rowMeans(mat_out_boot, na.rm = TRUE); rm(mat_out_boot)
        
        # store results to master file
        mat_out_tmp               <- data.frame(N1 = N1_grid)
        mat_out_tmp$N0            <- rep(N0, N1_grid_l)
        mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
        mat_out_tmp$name          <- "upstrap_power"
        mat_out_tmp$effect0       <- coef_x1
        mat_out_tmp$effect1       <- effect1
        mat_out_tmp$value         <- value
        mat_out_tmp$elapsed_s     <- NA
        
        mat_out_all <- rbind(mat_out_all, mat_out_tmp)
        rm(mat_out_tmp, value, fit)
    }
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH simr
    print("--- ESTIMATE POWER WITH simr")
    
    t_i_1 <- Sys.time()
    fit   <- lm(y ~ x1 + x2 + x3, data = datN0)
    # The along argument specifies which variable is being extended, 
    # and n specifies how many levels to replace it with
    fit_ext  <- simr::extend(fit, along = "subjid", n = N1_max)
    pc_out   <- simr::powerCurve(fit_ext, along = "subjid", breaks = N1_grid, nsim = B_boot, progress = FALSE)
    pc_out_s <- summary(pc_out)
    value    <- pc_out_s$mean
    t_i_2 <- Sys.time()
    t_i_12_diff <- round(as.numeric(difftime(t_i_2, t_i_1, units = "secs")), 1)
    # store results to master file
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$name          <- "simr_power"
    mat_out_tmp$effect0       <- coef_x1
    mat_out_tmp$effect1       <- NA
    mat_out_tmp$value         <- value
    mat_out_tmp$elapsed_s     <- t_i_12_diff
    mat_out_all <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value, fit, fit_ext, pc_out, pc_out_s, t_i_12_diff, t_i_1, t_i_2)
    
    
    # ------------------------------------------------------------------------------
    # ESTIMATE POWER WITH simr, modify an effect size 
    print("--- ESTIMATE POWER WITH simr, modify an effect size ")
    
    for (effect1 in effect1_grid){ # effect1 <- 0.3
        fit   <- lm(y ~ x1 + x2 + x3, data = datN0)
        coef(fit)["x1"] <- effect1
        fit_ext  <- simr::extend(fit, along = "subjid", n = N1_max)
        pc_out   <- simr::powerCurve(fit_ext, along = "subjid", breaks = N1_grid, nsim = B_boot, progress = FALSE)
        pc_out_s <- summary(pc_out)
        value    <- pc_out_s$mean
        # store results to master file
        mat_out_tmp               <- data.frame(N1 = N1_grid)
        mat_out_tmp$N0            <- rep(N0, N1_grid_l)
        mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
        mat_out_tmp$name          <- "simr_power"
        mat_out_tmp$effect0       <- coef_x1
        mat_out_tmp$effect1       <- effect1
        mat_out_tmp$value         <- value
        mat_out_tmp$elapsed_s     <- NA
        mat_out_all <- rbind(mat_out_all, mat_out_tmp)
        rm(mat_out_tmp, value, fit, fit_ext, pc_out, pc_out_s)
    }
    
    
    # ------------------------------------------------------------------------------
    # GET SINGLE RESULT
    print("--- GET SINGLE RESULT")
    
    mat_out_boot <- matrix(NA, nrow = N1_grid_l, ncol = 1)
    for (N1_grid_idx in 1 : length(N1_grid)){ # N1_grid_idx <-1
        N1 <- N1_grid[N1_grid_idx]
        datN0 <-  dat %>% filter(subjid_arm <= N1)
        fit  <- lm(y ~ x1 + x2 + x3, data = datN0)
        fit_s <- summary(fit)
        fit_coefpval  <- fit_s$coefficients[2, 4]
        mat_out_boot[N1_grid_idx, 1] <- (fit_coefpval < 0.05) * 1
    }
    value <- rowMeans(mat_out_boot, na.rm = TRUE)
    # store results to master file
    mat_out_tmp               <- data.frame(N1 = N1_grid)
    mat_out_tmp$N0            <- rep(N0, N1_grid_l)
    mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
    mat_out_tmp$name          <- "run_result"
    mat_out_tmp$effect0       <- coef_x1
    mat_out_tmp$effect1       <- NA
    mat_out_tmp$value         <- value
    mat_out_tmp$elapsed_s     <- NA
    mat_out_all <- rbind(mat_out_all, mat_out_tmp)
    rm(mat_out_tmp, value)
    
    
    # ------------------------------------------------------------------------------
    # GET SINGLE RESULT, modify an effect size 
    print("--- GET SINGLE RESULT, modify an effect size ")
    
    for (effect1 in effect1_grid){ # effect1 <- 0.3
        y_i_modified   <- coef_x0 + (effect1 * x1_i) + (coef_x2 * x2_i) + (coef_x3 * x3_i) + eps_i
        dat_modified   <- data.frame(y = y_i_modified, 
                                     x1 = x1_i, x2 = x2_i, x3 = x3_i,
                                     subjid = subjid_i, subjid_arm = subjid_arm_i)
        mat_out_boot <- matrix(NA, nrow = N1_grid_l, ncol = 1)
        for (N1_grid_idx in 1 : length(N1_grid)){ # N1_grid_idx <-1
            N1 <- N1_grid[N1_grid_idx]
            datN0 <-  dat_modified %>% filter(subjid_arm <= N1)
            fit  <- lm(y ~ x1 + x2 + x3, data = datN0)
            fit_s <- summary(fit)
            fit_coefpval  <- fit_s$coefficients[2, 4]
            mat_out_boot[N1_grid_idx, 1] <- (fit_coefpval < 0.05) * 1
        }
        value <- rowMeans(mat_out_boot, na.rm = TRUE)
        # store results to master file
        mat_out_tmp               <- data.frame(N1 = N1_grid)
        mat_out_tmp$N0            <- rep(N0, N1_grid_l)
        mat_out_tmp$arrayjob_idx  <- rep(arrayjob_idx, N1_grid_l)
        mat_out_tmp$name          <- "run_result"
        mat_out_tmp$effect0       <- coef_x1
        mat_out_tmp$effect1       <- effect1
        mat_out_tmp$value         <- value
        mat_out_tmp$elapsed_s     <- NA
        mat_out_all <- rbind(mat_out_all, mat_out_tmp)
        rm(mat_out_tmp, value)
    }
    
    # ------------------------------------------------------------------------------
    # SAVE TO FILE 
    out_fpath_raw <- paste0(res_fdir_raw, "/arrayjob_", arrayjob_idx, ".rds")
    saveRDS(object = mat_out_all, file = out_fpath_raw)
}


t2 <- Sys.time()
t2-t1

# number of repetitions of experiment 
# R_rep   <- 30
# number of boot repetitions within one experiment, one setup
# B_boot  <- 100
# Time difference of 12.25469 mins





