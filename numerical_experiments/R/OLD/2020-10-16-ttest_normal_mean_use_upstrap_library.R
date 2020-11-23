
rm(list = ls())

# libraries
library(data.table)
library(dplyr)
# devtools::install_github("ccrainic/upstrap")
# devtools::install_github("martakarass/upstrap")
library(upstrap)

# parameters definition
n0 <- 200
mu  <- 0.1
B_boot  <- 1000


## -----------------------------------------------------------------------------
## upstrap example

fname = system.file("extdata", "shhs1.txt", package = "upstrap")
dat = read.table(file = fname,  header = TRUE, na.strings="NA")

plot(dat$bmi_s1,dat$rdi4p,pch=".",col=rgb(0,0,1,.2),cex=3,ylim=c(0,50),bty="l",
     cex.axis=1.2,col.axis="blue",main=NULL,xlab="Body Mass Index",
     ylab="Respiratory Disturbance Index (4%)")
points(mean(dat$bmi_s1,na.rm=TRUE),mean(dat$rdi4p,na.rm=TRUE),cex=2,pch=19,col=rgb(1,0,0,0.5))

# calculate the mean BMI and RDI
round(mean(dat$bmi_s1,na.rm=TRUE),digits=2)
round(mean(dat$rdi4p,na.rm=TRUE),digits=2)

## Fit a basic regression model
# define the moderate to severe sleep apnea variable from rdi4p
MtS_SA=dat$rdi4p>=15

# add the MtS_SA variable to the data frame
dat$MtS_SA=MtS_SA
mean(MtS_SA)

fit<-glm(MtS_SA~gender+age_s1+bmi_s1+HTNDerv_s1+age_s1*HTNDerv_s1,
         family="binomial",data=dat)
summary(fit)



# set the seed for reproducibility
set.seed(08132018)

# sample size of the original data
n_oss=dim(dat)[1]

# number of grid points for the multiplier of the sample size
n_grid_multiplier=21

# minimum and maximum multiplier for the sample size
# here they are set to 1 (same sample size) and 5 (5 times the original sample size)
min_multiplier=1
max_multiplier=5

# set the grid of multipliers for the original sample size
# here 1.2 stands for a sampel size that is 20% larger than the original sample size
multiplier_grid=seq(from=min_multiplier, to=max_multiplier,length=n_grid_multiplier)

# vector of the new sample sizes
new_sample_sizes = n_oss * multiplier_grid

n_upstrap=100

statistic = function(data) {
  fit <- glm(MtS_SA~gender+age_s1+bmi_s1+HTNDerv_s1+age_s1*HTNDerv_s1, 
             family="binomial", 
             data=data)
  # could use broom
  #  tid = broom::tidy(fit)
  #  pval = tid$p.value[ tid$term == "HTNDerv_s1"]
  smod <- coef(summary(fit))
  pval <- smod["HTNDerv_s1", "Pr(>|z|)"]
  pval < 0.05
}

# invoke upstap for each new sample size
check = upstrap(dat, statistic = statistic, 
                R = n_upstrap, 
                new_sample_size = new_sample_sizes)

power_check = colMeans(check)
power_check

plot(multiplier_grid,power_check,type="l",col="blue",lwd=3,
     xlab="Factor by which the sample size is multiplied",
     ylab="Power to detect the HTN effect",
     bty="l",cex.axis=1.5,cex.lab=1.4,col.lab="blue",
     col.axis="blue",main=NULL)
# add horizontal lines to indicate power equal to 0.8 (orange) and 0.9 (red).
abline(h=0.8,lwd=3,col="orange")
abline(h=0.9,lwd=3,col="red")



# -----------------------------------------------------

rm(list = ls())

library(data.table)
library(dplyr)
library(upstrap)

# parameters definition
n0        <- 200
mu        <- 0.1
R_boot    <- 1000
power_val <- 0.8

set.seed(1)

sample_obs <- rnorm(n = n0, mean = mu, sd = 1)

## Approach 1: power.t.test
out   <- power.t.test(delta = mean(sample_obs), sd = sd(sample_obs), power = 0.8)
out_n <- ceiling(out$n)
out_n

## Approach 2: upstrap
rm(statistic, new_sample_sizes)
statistic = function(data) {
  fit <- t.test(data)
  pval <- fit$p.value
  pval < 0.05
}
new_sample_size <- round(seq(n0, round(out_n * 1.5), length.out = 20))
data <- sample_obs

upstrap_out <- upstrap(
  data = sample_obs,
  statistic = statistic,
  R = R_boot,
  new_sample_size = new_sample_size)

colMeans(upstrap_out)
