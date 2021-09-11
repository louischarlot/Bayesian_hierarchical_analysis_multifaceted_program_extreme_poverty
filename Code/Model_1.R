# Author: Louis Charlot
# Title: Model 1 (See master thesis for more details)

# This is our code for Model 1. 
# It begins with a Setup that downloads the needed libraries.
# Then, it downloads the data of the multifaceted program and prepares it.
# Then, it executes the site-level regressions and saves the estimates and associated standard errors.
# Finally it calls the hierarchical model from Model_1.stan and displays the results.





#####################################################################################################################
##### SETUP #########################################################################################################
#####################################################################################################################
# We downloads the needed libraries.

# Set working directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### Preliminaries ###################################################################################################

rm(list = ls())

chooseCRANmirror(ind=90)


# Installs packages if not already installed, then loads packages 
installer <- TRUE
if(installer == TRUE){
  install.packages('foreign')
  install.packages('Hmisc')
  install.packages('xtable')
  install.packages('coda')
  install.packages('stargazer')
  install.packages('dummies')
  install.packages('zoo')
  install.packages('lmtest')
  install.packages('plm')
  install.packages('sandwich')
  install.packages('Matrix')
  install.packages("ggplot2")
  install.packages("gridExtra")
  install.packages("shinystan")
  install.packages("tidyverse")
  install.packages("readr")
  install.packages("foreign")
  install.packages("readstata13")    # read original Stata file
  install.packages("loo")  # compare Stan models
  install.packages("parallel")  # to detect number of cores
  install.packages("dplyr")
  install.packages("bayesplot")
}

library(dummies)
library(multiwayvcov)
library(stargazer)
library(foreign)
library(Hmisc)
library(xtable)
library(coda)
library(zoo)
library(plm)
library(lmtest)
library(sandwich)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(rstan)
library(shinystan)
library(tidyverse)
library(readr)
library(foreign)
library(readstata13)
library(loo)
library(parallel)
library(dplyr)
library(bayesplot)
library(ggplot2)



#####################################################################################################################
##### CHOICE OF VARIABLE OF INTEREST ################################################################################
#####################################################################################################################
y_interest <- 'asset_index_end'
#####################################################################################################################
#####################################################################################################################


## load an prepare data ################################################################################################


# We load the data:
# The "version" data_1 only converts "Yes" into "1" compared to data_2.
data_1 <- read.dta13("pooled_hh_postanalysis.dta")
data_2 <- read.dta13("pooled_hh_postanalysis.dta", nonint.factors=T,generate.factors=T)


# We create a function that removes rows for which a specified column has NAs:
Remove_rows_with_NA <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Stan cannot handle missing values in data automatically, so no element of the data we use can contain NA values.
# We remove rows where the variables we will use are NAs:
data_1 <- Remove_rows_with_NA(data_1, c("treatment",          
                                        "asset_index_end"))
data_2 <- Remove_rows_with_NA(data_2, c("treatment",          
                                        "asset_index_end"))


# Data by country:
data_Ethiopia  <- data_2[data_2$country=='Ethiopia',]
data_Ghana  <- data_2[data_2$country=='Ghana',]
data_Honduras  <- data_2[data_2$country=='Honduras',]
data_India  <- data_2[data_2$country=='India (Bandhan)',]
data_Pakistan  <- data_2[data_2$country=='Pakistan',]
data_Peru  <- data_2[data_2$country=='Peru',]


#####################################################################################################################
##### SOME INTRODUCTORY ANALYSIS ####################################################################################
#####################################################################################################################

regression_equation <- asset_index_end ~ treatment


# Site-level regressions for each site s:
regression_comsumption_data_Ethiopia <- summary(lm(regression_equation , data = data_Ethiopia))
print(regression_comsumption_data_Ethiopia)
regression_comsumption_data_Ghana <- summary(lm(regression_equation , data = data_Ghana))
print(regression_comsumption_data_Ghana)
regression_comsumption_data_Honduras <- summary(lm(regression_equation , data = data_Honduras))
print(regression_comsumption_data_Honduras)
regression_comsumption_data_India <- summary(lm(regression_equation , data = data_India))
print(regression_comsumption_data_India)
regression_comsumption_data_Pakistan <- summary(lm(regression_equation , data = data_Pakistan))
print(regression_comsumption_data_Pakistan)
regression_comsumption_data_Peru <- summary(lm(regression_equation , data = data_Peru))
print(regression_comsumption_data_Peru)



# Save the estimate hat(tau_s) and its associated standard error hat(sigma_s) for each regression (= for each site) :
hat_tau_s_Ethiopia <- regression_comsumption_data_Ethiopia$coefficients[2,1]
sigma_s_Ethiopia <- regression_comsumption_data_Ethiopia$coefficients[2,2]
hat_tau_s_Ghana <- regression_comsumption_data_Ghana$coefficients[2,1]
sigma_s_Ghana <- regression_comsumption_data_Ghana$coefficients[2,2]
hat_tau_s_Honduras <- regression_comsumption_data_Honduras$coefficients[2,1]
sigma_s_Honduras <- regression_comsumption_data_Honduras$coefficients[2,2]
hat_tau_s_India <- regression_comsumption_data_India$coefficients[2,1]
sigma_s_India <- regression_comsumption_data_India$coefficients[2,2]
hat_tau_s_Pakistan <- regression_comsumption_data_Pakistan$coefficients[2,1]
sigma_s_Pakistan <- regression_comsumption_data_Pakistan$coefficients[2,2]
hat_tau_s_Peru <- regression_comsumption_data_Peru$coefficients[2,1]
sigma_s_Peru <- regression_comsumption_data_Peru$coefficients[2,2]

hat_tau_s <- c(hat_tau_s_Ethiopia, hat_tau_s_Ghana, hat_tau_s_Honduras, hat_tau_s_India, hat_tau_s_Pakistan, hat_tau_s_Peru)
sigma_s <- c(sigma_s_Ethiopia, sigma_s_Ghana, sigma_s_Honduras, sigma_s_India, sigma_s_Pakistan, sigma_s_Peru)








#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
##### BAYESIAN HIERARCHICAL ANALYSIS (STAN) #########################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# We now call the hierarchical model from Model_1.stan and display the results.
# Look Master Thesis for more details.

# Preparation of the data for Stan:
data_for_stan_1 <- list(
  S = 6, # number of sites
  hat_tau_s = hat_tau_s,
  sigma_s = sigma_s
)

# We call the stan function to draw posterior samples:
fit_1_optimized <- stan(
  file = "Model_1.stan",  # Stan model
  data = data_for_stan_1, # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = max(2000/10, 1),  # default actualization of progress shown
  control=list(adapt_delta=0.99) # We reduce the step size.
)

plot(fit_1_optimized, pars = c("tau", "sigma", "tau_s"))


# To lanch shinystan:
#launch_shinystan(fit_1_optimized)



# BAYESPLOT: #####################################################################################################################
# To have nice graphical representations.

color_scheme_set("brightblue")
posterior_1_optimized <- as.matrix(fit_1_optimized)
plot_title <- ggtitle("",expression(paste("Model 1 with priors ", tau %~% N(0,5), " and ", sigma %~% halfCauchy(0,5))))

# Density plots computed from posterior draws with all chains merged, with uncertainty intervals shown as shaded areas under the curves.
plot_post_distrib_1_optimized <- mcmc_areas(posterior_1_optimized,
                                                      pars = c("tau_s[1]", "tau_s[2]", "tau_s[3]", "tau_s[4]", "tau_s[5]", "tau_s[6]"),
                                                      point_est = c("none"),
                                                      prob = 0.95) + 
  plot_title +
  scale_y_discrete(limits=c("tau_s[6]", "tau_s[5]", "tau_s[4]", "tau_s[3]", "tau_s[2]", "tau_s[1]"),
                   labels=c(    # Labels fot the plot
                     expression(paste(tau[6], " (Peru ITT)")),
                     expression(paste(tau[5], " (Pakistan ITT)")),
                     expression(paste(tau[4], " (India ITT)")), 
                     expression(paste(tau[3], " (Honduras ITT)")),
                     expression(paste(tau[2], " (Ghana ITT)")), 
                     expression(paste(tau[1], "  (Ethiopia ITT)"))
                   )) + 
  xlim(-0.1, 1)
ggplot2::ggsave("plot_Model_1_optimized.png", width = 6, height = 6)
plot_post_distrib_1_optimized


# And for sigma:
plot_sigma_post_distrib_1_optimized <- mcmc_areas(posterior_1_optimized,
                                            pars = c("sigma"),
                                            point_est = c("none"),
                                            prob = 0.95) + 
  plot_title +
  xlim(-0.1, 1)
ggplot2::ggsave("plot_sigma_Model_1_optimized.png", width = 6, height = 4)
plot_sigma_post_distrib_1_optimized








