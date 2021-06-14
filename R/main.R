# Bayesian time-varying autoregressive models of COVID-19 epidemics

rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
setwd("/Users/barbaratarantino/Desktop/BAYES_DynPoisson")
library(fda)
library(pracma)
library(scoringutils)
library(covidHubUtils)
library(bayestestR)
library(lubridate)
library(dplyr)
library(tscount)
library(tidyverse)
source("HMC.R")
source("BAYES_log.R")
source("tvBAYES_log.R")
source("posterior_help.R")
source("main_functions.R")

# Select country of interest and import data -----------------------------------

df <- read.data("df_bayes_all.csv", "Italy")

# Train and test (last 7 days) -------------------------------------------------

data_list <- split(df, npi_lag = NULL)

# Fit the model ----------------------------------------------------------------

model = "ARX"
data_train = data_list[["train"]][["data_train"]]
xreg_train = data_list[["train"]][["xreg_train"]]

fit <- fit.ARX.log(data =data_train, Xreg = xreg_train, order=1,
                    Total_itr = 10000, burn=5000)

saveRDS(fit, "fit_ARX_log.RDS")

# Posterior diagnostics --------------------------------------------------------

fit_post <- bayes.posterior(fit, data_list, 
                            model = list(past_obs = 1, past_mean = NULL, tv = FALSE), 
                            link = "log", distr = "poisson", xreg = TRUE, 
                            country = as.character(unique(df$Country)))

fit_pred <- bayes.predict(fit, data_list, 
                          model = list(past_obs = 1, past_mean = NULL, tv = FALSE), 
                          coeff_list = fit_post[["coefficients"]][["coefficient.distribution"]], 
                          n.ahead = 7, validation = TRUE, xreg = TRUE, 
                          country = as.character(unique(df$Country)))

