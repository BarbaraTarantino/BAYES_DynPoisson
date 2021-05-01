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
  dplyr::select(Number) %>%
  pull()
xreg_train <- df_train %>%
  dplyr::select(starts_with(c("c", "h"))) %>%
  dplyr::select(!Country)

df_test = df %>% tail(7)

# Fit the model

fit <- fit.INGARCHX.log(data_train,xreg_train)

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
                   data= data, order_list = list(order1=1, order2=1))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At)
  
  vt <- predict(data=data, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data, date = df %>% pull(Date), 
                          coeff_list=coeff_list, 
                          order_list = list(order1 = 1, order2 = 0),
                          n.ahead=1)
  
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
                   data= data, order_list = list(order1=1, order2=0))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data, order_list = list(order1 = 1, order2 = 0))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data, xreg=Xreg, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  vt <- predict(data=data, xreg = Xreg, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data, date = df %>% pull(Date), 
                         xreg = Xreg, coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 0),
                         n.ahead=1)
  

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
      
      bsplines <- bsplines_fun(data, order_list = list(order1 = 1, order2 = 0),
                               knot = 4, norder = 4)
      
      mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                       data= data, order_list = list(order1=1, order2=0),
                       bsplines = bsplines)
      
      At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                        data = data, order_list = list(order1 = 1, order2 = 0), 
                        bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
      
      coeff_list <- list(mut=mut, At=At)
      
      vt <- predict(data=data, coeff_list=coeff_list, 
                    order_list = list(order1 = 1, order2 = 0))
      
      fcst <- forecast_ahead(data=data, date = df %>% pull(Date), 
                              coeff_list=coeff_list, 
                              order_list = list(order1 = 1, order2 = 0),
                              n.ahead=1)
      
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
  
  bsplines <- bsplines_fun(data, order_list = list(order1 = 1, order2 = 0),
                           knot = 4, norder = 4)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data, order_list = list(order1=1, order2=0),
                   bsplines = bsplines)
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data, order_list = list(order1 = 1, order2 = 0), 
                    bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data, xreg=Xreg, order_list = list(order1 = 1, order2 = 0),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  vt <- predict(data=data, xreg = Xreg, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 0))
  
  fcst <- forecast_ahead(data=data, date = df %>% pull(Date), 
                          xreg = Xreg, coeff_list=coeff_list, 
                          order_list = list(order1 = 1, order2 = 0),
                          n.ahead=1)
  
} else if (model == "INGARCH"){
  
} else if (model == "INGARCHX"){
  
} else if (model == "tvINGARCH"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "sigma2latfn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaKfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  bsplines <- bsplines_fun(data, order_list = list(order1 = 1, order2 = 1),
                           knot = 4, norder = 4)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data, order_list = list(order1=1, order2=1),
                   bsplines_list = bsplines)
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data, order_list = list(order1 = 1, order2 = 1), 
                       bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]])
  
  vt <- predict(data=data, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data, date = df %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=1)
  
} else if (model == "tvINGARCHX"){
  
  # traceplot(fit, "deltamufn", tv = "TRUE")
  # traceplot(fit, "sigma2latfn", tv = "TRUE")
  # traceplot(fit, "deltaRfn", tv = "TRUE")
  # traceplot(fit, "deltaKfn", tv = "TRUE")
  # traceplot(fit, "deltaCfn", tv = "TRUE")
  # traceplot(fit, "pred", tv = "TRUE")
  # traceplot(fit, "vtp", tv = "TRUE")
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  bsplines <- bsplines_fun(data)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data, order = list(order1=1, order2=1),
                   bsplines = bsplines)
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data, order = list(order1 = 1, order2 = 1), 
                       bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data, xreg=Xreg, order = list(order1 = 1, order2 = 1),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]], Ct=Ct)
  
  # Predictions 
  
  vt <- predict(data=data, xreg=Xreg, coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  fcst <- forecast_ahead(data=data, xreg =Xreg, 
                         date = df %>% pull(Date), coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=1)
  
}
  
  
  
  
