
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Upstrap manuscript

This GitHub directory contains code for all simulations and real data
application example presented in the manuscript draft “Upstrap for
Estimating Power and Sample Size in Complex Models”.

## Simulations

The Joint High Performance Computing Exchange
([JHPCE](https://jhpce.jhu.edu/)) is a High-Performance Computing (HPC)
facility in the Department of Biostatistics at the Johns Hopkins
Bloomberg School of Public Health was used to execute and aggregate the
numerical experiments.

The aggregated data output files are stored within this repository and
can be used to reproduce tables and files included in the manuscript.

### R code scripts ran on the JHPCE cluster

-   [run\_onesample\_ttest.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_onesample_ttest.R)
    – R script to run code for simulation problem 1
-   [run\_twosample\_ttest.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_twosample_ttest.R)
    – R script to run code for simulation problem 2
-   [run\_lm\_testcoef.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lm_testcoef.R)
    – R script to run code for simulation problem 3
-   [run\_glm\_testcoef.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_glm_testcoef.R)
    – R script to run code for simulation problem 4
-   [run\_lmm\_testcoef.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lmm_testcoef.R)
    – R script to run code for simulation problem 5
-   [run\_glmm\_testcoef.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_glmm_testcoef.R)
    – R script to run code for simulation problem 6 </br></br>
-   [run\_twosample\_ttest\_vary\_covprop.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_twosample_ttest_vary_covprop.R)
    – R script to run code for simulation problem 7
-   [run\_lm\_testcoef\_vary\_covprop.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lm_testcoef_vary_covprop.R)
    – R script to run code for simulation problem 8
-   [run\_lmm\_testcoef\_vary\_covprop.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lmm_testcoef_vary_covprop.R)
    – R script to run code for simulation problem 9 </br></br>
-   [run\_twosample\_ttest\_vary\_nobs\_sd.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_twosample_ttest_vary_nobs_sd.R)
    – R script to run code for simulation problem 10
-   [run\_lm\_testcoef\_vary\_nobs\_sd.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lm_testcoef_vary_nobs_sd.R)
    – R script to run code for simulation problem 11
-   [run\_lmm\_testcoef\_vary\_nobs\_sd.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/run_lmm_testcoef_vary_nobs_sd.R)
    – R script to run code for simulation problem 12</br></br>
-   [agg\_sim\_123456.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/agg_sim_123456.R)
    – R script to aggregate results from independent experiment
    repetitions for simulation problems 1-6
-   [agg\_sim\_vary\_covprop.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/agg_sim_vary_covprop.R)
    – R script to aggregate results from independent experiment
    repetitions for simulation problems 7-9
-   [agg\_sim\_vary\_nobs\_sd.R](https://github.com/martakarass/upstrap_manuscript/blob/master/numerical_experiments/R/agg_sim_vary_nobs_sd.R)
    – R script to aggregate results from independent experiment
    repetitions for simulation problems 10-12

### Aggregated simulation results

-   [/results\_CL\_shared](https://github.com/martakarass/upstrap_manuscript/tree/master/numerical_experiments/results_CL_shared)
    – directory with files containing aggregated results from
    independent experiment repetitions for simulation problems 1-12

### R code scripts – figures

-   [/results\_CL\_shared](https://github.com/martakarass/upstrap_manuscript/tree/master/numerical_experiments/results_CL_shared)
    – directory with files containing aggregated results from
    independent experiment repetitions for simulation problems 1-12

### R code scripts – tables
