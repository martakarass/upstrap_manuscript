#' @description 
#' Script to download from the cluster results of experiment in which we 
#' n_rep = 100 x 1,000 times simulated data from normal distribution and 
#' computed sample size needed to observe effect being statistically 
#' significant, using: 
#' (a) t.test() approach
#' (b) upstrap approach 
#' 
#' @author 
#' Marta Karas <mkaras2@jhu.edu>


## -----------------------------------------------------------------------------

# ssh -X mkaras@jhpce01.jhsph.edu
# cd /users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean
# qrsh -l mem_free=30G,h_vmem=30G,h_stack=256M
# module load conda_R
# R

library(data.table)
library(dplyr)
fdir <- "/users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean/"
list.files(fdir)

fnames_powerttest <- list.files(fdir)[grepl(pattern = "powerttest", list.files(fdir))]
fnames_upstrap    <- list.files(fdir)[grepl(pattern = "upstrap", list.files(fdir))]

t1 <- Sys.time()
out_powerttest_dt <- rbindlist(lapply(paste0(fdir, fnames_powerttest), fread))
t2 <- Sys.time()
t2-t1
# Time difference of 7.067129 secs

t1 <- Sys.time()
out_upstrap_dt <- rbindlist(lapply(paste0(fdir, fnames_upstrap), fread))
t2 <- Sys.time()
t2-t1
# Time difference of 6.500458 secs

# turn data to data frame 
out_powerttest <- as.data.frame(out_powerttest_dt)
out_upstrap    <- as.data.frame(out_upstrap_dt)

out_powerttest <- mutate(out_powerttest, fileid = floor((row_number()-1) / 1000), .before = everything())
out_upstrap    <- mutate(out_upstrap, fileid = floor((row_number()-1) / 1000), .before = everything())

out_powerttest <- mutate(out_powerttest, approach = "powerttest", .before = everything())
out_upstrap    <- mutate(out_upstrap, approach = "upstrap", .before = everything())

# combine two files together
out <- rbind(out_powerttest, out_upstrap); rm(out_powerttest, out_upstrap)

# aggregate across (approach, file_id) 
out_agg_across  <- 
  out %>%  
  group_by(approach, file_id) %>% 
  summarise(across(.cols = everything(), 
                   .fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               median = ~median(.x, na.rm = TRUE))))
dim(out_agg_across)

# aggregate across (approach) 
out_agg_all  <- 
  out %>%  
  select(-file_id) %>%
  group_by(approach) %>% 
  summarise(across(.cols = everything(), 
                   .fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               median = ~median(.x, na.rm = TRUE))))
dim(out_agg_all)

# save to file 
saveRDS(out_agg_across, paste0(fdir, "out_agg_across"))
saveRDS(out_agg_all, paste0(fdir, "out_agg_all"))

# sftp mkaras@jhpce-transfer01.jhsph.edu
# get /users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean/out_agg_across /Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean/out_agg_across
# get /users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean/out_agg_all /Users/martakaras/Dropbox/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL/2020-11-23-ttest_normal_mean/out_agg_all







