rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
setwd("/Users/barbaratarantino/Desktop/BAYES_DynPoisson")
library(fda)
library(pracma)
library(dplyr)
library(tscount)
library(tidyverse)
source("HMC.R")
source("BAYES_log.R")
source("tvBAYES_log.R")
source("posterior_help.R")

# Select country of interest
# Import count data + covariates

c = "Italy"
df <- read.csv("df_bayes.csv") %>%
  filter(Country == c)

# Train and test (last 7 days)

df_train = df %>% slice(1:(n()-7))
data_train <- df_train %>% 
  mutate(Number = as.numeric(Number)) %>%
  dplyr::select(Number) %>%
  pull()
xreg_train <- df_train %>%
  dplyr::select(starts_with(c("c", "h"))) %>%
  dplyr::select(!Country) %>%
  mutate_all(as.numeric)

df_test = df %>% tail(7)

# Fit the model

fit1 <- fit.ARX.log(data =data_train, Xreg = xreg_train, order=1,
                   Total_itr = 10000, burn=5000)
fit2 <- fit.INGARCHX.log(data_train, Xreg = xreg_train, order1=1, order2 =1,
                        Total_itr = 10000, burn=5000)
fit3 <- fit.tvARX.log(data=data_train, Xreg = xreg_train, knot = 4, norder=4,
                     order = 1, Total_itr = 10000, burn=5000)
fit4 <- fit.tvINGARCHX.log(data_train, Xreg = xreg_train, knot = 4, norder=4,
                          order1 = 1, order2 = 1, Total_itr = 10000, burn=5000)

saveRDS(fit1, "fit_ARX_log.RDS")
saveRDS(fit2, "fit_INGARCHX_log.RDS")
saveRDS(fit3, "fit_tvARX_log.RDS")
saveRDS(fit4, "fit_tvINGARCHX_log.RDS")

# Posterior diagnostics 

model = "tvAR" 

if (model == "AR"){
  
  # traceplot(fit, "deltamufn", tv = "FALSE")
  # traceplot(fit, "deltaRfn", tv = "FALSE")
  # traceplot(fit, "pred", tv = "FALSE")
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=1))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At)
  
  vt <- predict(data=data_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                          coeff_list=coeff_list, 
                          order_list = list(order1 = 1, order2 = 0),
                          n.ahead=7)
  
} else if (model == "ARX"){

  # traceplot(fit, "deltamufn", tv = "FALSE")
  # traceplot(fit, "deltaRfn", tv = "FALSE")
  # traceplot(fit, "deltaCfn", tv = "FALSE")
  # traceplot(fit, "pred", tv = "FALSE")
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=0))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xxreg_train, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  vt <- predict(data=data_train, xreg = xreg_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         xreg = xreg_train, coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 0),
                         n.ahead=7)
  

}else if (model == "tvAR"){
    
    # traceplot(fit, "deltamufn", tv = "TRUE")
    # traceplot(fit, "deltaRfn", tv = "TRUE")
    # traceplot(fit, "pred", tv = "TRUE")
    # traceplot(fit, "vtp", tv = "TRUE")
      
      delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
      p_summary = list()
      for (j in delta){
        sum <- summary(fit, j)
        p_summary[[j]] <- sum
      }
      
      bsplines <- bsplines_fun(data_train, order_list = list(order1 = 1, order2 = 0),
                               knot = 4, norder = 4)
      
      mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                       data= data_train, order_list = list(order1=1, order2=0),
                       bsplines_list =list(bsplines=bsplines, knot = 4, norder=4))
      
      At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                        data = data_train, order_list = list(order1 = 1, order2 = 0), 
                        bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
      
      coeff_list <- list(mut=mut, At=At)
      
      vt <- predict(data=data_train, coeff_list=coeff_list, 
                    order_list = list(order1 = 1, order2 = 0))
      
      fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                              coeff_list=coeff_list, 
                              order_list = list(order1 = 1, order2 = 0),
                              n.ahead=7)
      
} else if(model == "tvARX"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaCfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  bsplines <- bsplines_fun(data_train, order_list = list(order1 = 1, order2 = 0),
                           knot = 4, norder = 4)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=0),
                   bsplines = bsplines)
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0), 
                    bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order_list = list(order1 = 1, order2 = 0),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  vt <- predict(data=data_train, xreg = xreg_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                          xreg = xreg_train, coeff_list=coeff_list, 
                          order_list = list(order1 = 1, order2 = 0),
                          n.ahead=7)
  
} else if (model == "INGARCH"){
  
  # traceplot(fit, "deltamufn", tv = "FALSE")
  # traceplot(fit, "sigma2latfn", tv = "FALSE")
  # traceplot(fit, "deltaRfn", tv = "FALSE")
  # traceplot(fit, "deltaKfn", tv = "FALSE")
  # traceplot(fit, "pred", tv = "FALSE")
  
  for (i in 1:length(fit[["deltaRfn"]])){
    deltaRfn <- unlist(fit[["deltaRfn"]][[i]])
    deltaKfn <- unlist(fit[["deltaKfn"]][[i]])
    fit[["deltaRfn"]][[i]] <- (deltaRfn+deltaKfn)/2
    fit[["deltaKfn"]][[i]] <- (deltaRfn-deltaKfn)/2
  }
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=1))
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data_train, order_list = list(order1 = 1, order2 = 1))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]])
  
  vt <- predict(data=data_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
} else if (model == "INGARCHX"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "sigma2latfn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaKfn", tv = "TRUE")
  # traceplot(fit, "deltaCfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  # traceplot(fit, "vtp", tv = "TRUE")
  
  for (i in 1:length(fit[["deltaRfn"]])){
    deltaRfn <- unlist(fit[["deltaRfn"]][[i]])
    deltaKfn <- unlist(fit[["deltaKfn"]][[i]])
    fit[["deltaRfn"]][[i]] <- (deltaRfn+deltaKfn)/2
    fit[["deltaKfn"]][[i]] <- (deltaRfn-deltaKfn)/2
  }
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order = list(order1=1, order2=1))
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data_train, order = list(order1 = 1, order2 = 1))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order = list(order1 = 1, order2 = 1))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]], Ct=Ct)
  
  # Predictions 
  
  vt <- predict(data=data_train, xreg=xreg_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data_train, xreg =xreg_train, 
                         date = df_train %>% pull(Date), coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
} else if (model == "tvINGARCH"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "sigma2latfn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaKfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  
  for (i in 1:length(fit[["deltaRfn"]])){
    deltaRfn <- unlist(fit[["deltaRfn"]][[i]])
    deltaKfn <- unlist(fit[["deltaKfn"]][[i]])
    fit[["deltaRfn"]][[i]] <- (deltaRfn+deltaKfn)/2
    fit[["deltaKfn"]][[i]] <- (deltaRfn-deltaKfn)/2
  }
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  bsplines <- bsplines_fun(data_train, order_list = list(order1 = 1, order2 = 1),
                           knot = 4, norder = 4)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=1),
                   bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data_train, order_list = list(order1 = 1, order2 = 1), 
                       bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]])
  
  vt <- predict(data=data_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
} else if (model == "tvINGARCHX"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "sigma2latfn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaKfn", tv = "TRUE")
  # traceplot(fit, "deltaCfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  # traceplot(fit, "vtp", tv = "TRUE")
  
  for (i in 1:length(fit[["deltaRfn"]])){
    deltaRfn <- unlist(fit[["deltaRfn"]][[i]])
    deltaKfn <- unlist(fit[["deltaKfn"]][[i]])
    fit[["deltaRfn"]][[i]] <- (deltaRfn+deltaKfn)/2
    fit[["deltaKfn"]][[i]] <- (deltaRfn-deltaKfn)/2
  }
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  bsplines <- bsplines_fun(data_train)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order = list(order1=1, order2=1),
                   bsplines = bsplines)
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data_train, order = list(order1 = 1, order2 = 1), 
                       bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order = list(order1 = 1, order2 = 1),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]], Ct=Ct)
  
  # Predictions 
  
  vt <- predict(data=data_train, xreg=xreg_train, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data_train, xreg =xreg_train, 
                         date = df_train %>% pull(Date), coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
}
  
  
  
  
