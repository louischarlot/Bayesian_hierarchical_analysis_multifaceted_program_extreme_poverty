# Author: Louis Charlot
# Title: Model 2bis (See master thesis for more details)

# This is our code for Model 2bis. 
# It begins with a Setup that downloads the needed libraries.
# Then, it downloads the data of the multifaceted program and prepares it.
# Then, it calls the hierarchical model from Model_2.stan and displays the results.



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
installer <- FALSE
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
# The "version" data_1 converts "Yes" into "1" compared to data_2.
data_1 <- read.dta13("pooled_hh_postanalysis.dta")
data_2 <- read.dta13("pooled_hh_postanalysis.dta", nonint.factors=T,generate.factors=T)

# We create a function that removes rows for which a specified column has NAs:
Remove_rows_with_NA <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Stan cannot handle missing values in data automatically, so no element of the data can contain NA values.
# We remove rows where the variables we will use are NAs:
data_1 <- Remove_rows_with_NA(data_1, c("treatment",          
                                        "asset_index_end", "asset_index_bsl"))

data_2 <- Remove_rows_with_NA(data_2, c("treatment",          
                                        "asset_index_end", "asset_index_bsl"))


# Data by country:
data_Ethiopia  <- data_2[data_2$country=='Ethiopia',]
data_Ghana  <- data_2[data_2$country=='Ghana',]
data_Honduras  <- data_2[data_2$country=='Honduras',]
data_India  <- data_2[data_2$country=='India (Bandhan)',]
data_Pakistan  <- data_2[data_2$country=='Pakistan',]
data_Peru  <- data_2[data_2$country=='Peru',]







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


#####################################################################################################################
##### 2bis) Second Model with baseline (STAN) #######################################################################
#####################################################################################################################
# We call the hierarchical model from Model_2.stan and display the results.
# Look Master Thesis for more details.



# Prepare the values that the sampler will use:
y <- data_1$asset_index_end # Vector of individual-level outcome of interest.
T <- data_1$treatment
N <- length(y)
S <- 6  # Number of sites.
site <- data_1$country #factor variable to split them out into K sites
# Individual-level predictors:
data_1$X_intercept <- rep(1,N)
# => 1 = Ethiopia, 2 = Ghana, 3 = Honduras, 4 = India, 5 = Pakistan, 6 = Peru
Covariates_X <- dplyr::select(data_1, X_intercept, treatment, asset_index_bsl) # also baseline values
X <- as.matrix(Covariates_X) # Matrix of individual-level predictors.
I <- length(Covariates_X) # Number of individual-level predictors.
# Site-level predictors:
Z_intercept <- rep(1,S)
site_health_component <- c(0, 1, 1, 1, 1, 1)
site_value_asset_transfer <- c(7.98, 6, 4.75, 6.53, 3.75, 17.14) # in local goat price
Covariates_Z <- data.frame(Z_intercept, site_health_component, site_value_asset_transfer)
Z <- as.matrix(Covariates_Z) # Matrix of site-level predictors.
tZ <- t(Z) # We transpose the group predictors.
J <- length(Covariates_Z)  # Number of site-level predictors.

# Data prepared for Stan:
data_for_stan <- list(N = N, 
                      S = S, 
                      I = I,
                      J = J,
                      y = y,       
                      X = X,
                      tZ = tZ,
                      site = site)



fit_Cholesky_Intercept_and_Baseline <- stan(
  file = "Model_2.stan",  # Stan program
  data = data_for_stan,    # data prepared for Stan
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = max(2000/10, 1),  # default actualization of progress shown
  control=list(adapt_delta=0.99, max_treedepth = 14) # We reduce step size.
)


print(fit_Cholesky_Intercept_and_Baseline, pars = c("beta[2,1]", "beta[2,2]", "beta[2,3]", "beta[2,4]", "beta[2,5]", "beta[2,6]"))

plot(fit_Cholesky_Intercept_and_Baseline,
     pars = c("beta[2,1]", "beta[2,2]", "beta[2,3]", "beta[2,4]", "beta[2,5]", "beta[2,6]"), outer_level = 0.95, ci_level = 0.5)




# BAYESPLOT: #####################################################################################################################
# To display beautifully the results:

color_scheme_set("brightblue")
posterior_Intercept_and_Baseline <- as.matrix(fit_Cholesky_Intercept_and_Baseline)
plot_title <- ggtitle("",expression(atop(paste("Model 2bis with priors ", gamma[k] %~% N(0,5), ", ", theta[k] %~% halfCauchy(0,2.5)),
                                         paste(" and ", Omega [L] %~% LKJcorrCholesky(2)))))


# Density plots computed from posterior draws with all chains merged, with uncertainty intervals shown as shaded areas under the curves.
plot_post_distrib_Intercept_and_Baseline <- mcmc_areas(posterior_Intercept_and_Baseline,
                                pars = c("beta[2,1]", "beta[2,2]", "beta[2,3]", "beta[2,4]", "beta[2,5]", "beta[2,6]"),
                                point_est = c("none"),
                                prob = 0.95) + 
  plot_title +
  scale_y_discrete(limits=c("beta[2,6]", "beta[2,5]", "beta[2,4]", "beta[2,3]", "beta[2,2]", "beta[2,1]"),
                   labels=c(    # Labels fot the plot
                     expression(paste(beta[2][6], " (Peru ITT)")),
                     expression(paste(beta[2][5], " (Pakistan ITT)")),
                     expression(paste(beta[2][4], " (India ITT)")), 
                     expression(paste(beta[2][3], " (Honduras ITT)")),
                     expression(paste(beta[2][2], " (Ghana ITT)")), 
                     expression(paste(beta[2][1], "  (Ethiopia ITT)"))
                   )) + 
  xlim(-0.1, 1)
ggplot2::ggsave("plot_Model_2bis_Cholesky_Intercept_and_Baseline.png", width = 6, height = 6)
plot_post_distrib_Intercept_and_Baseline






