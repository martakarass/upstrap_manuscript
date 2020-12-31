# ssh -X mkaras@jhpce01.jhsph.edu

# cd /users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/results_CL
# qrsh -l mem_free=10G,h_vmem=10G,h_stack=256M
# module load conda_R
# R

# cd /users/mkaras/_PROJECTS/upstrap_manuscript/numerical_experiments/bash_CL/2020-12-29-lmm_trt
# qsub -cwd -N ttest -t 1-100 ttest_normal_mean

library(here)
library(tidyverse)


## -----------------------------------------------------------------------------
## (3) LMM trt 

# out_agg_fname <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2020-12-29-lmm_trt")
out_agg_fname <- paste0(here::here(), "/numerical_experiments/results_CL_shared/2020-12-30-lmm_trt")
out_fdir <- paste0(here::here(), "/numerical_experiments/results_CL/2020-12-30-lmm_trt/")

out_fnames <- list.files(out_fdir, full.names = TRUE)
out_f <- do.call("rbind", lapply(out_fnames, readRDS))
dim(out_f)
str(out_f)

out_f_agg <-
  out_f %>% 
  group_by(N0, N1) %>% 
  mutate(cnt = n()) %>%
  group_by(N0, N1, cnt) %>% 
  summarise(across(.cols = starts_with("out_"), 
                   .fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               qt5 = ~quantile(.x, na.rm = TRUE, probs = 0.05),
                               qt25 = ~quantile(.x, na.rm = TRUE, probs = 0.25),
                               qt75 = ~quantile(.x, na.rm = TRUE, probs = 0.75),
                               qt95 = ~quantile(.x, na.rm = TRUE, probs = 0.95)))) %>%
  as.data.frame()
saveRDS(out_f_agg, out_agg_fname)