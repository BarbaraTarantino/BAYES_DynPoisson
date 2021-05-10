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
source("BAYES_log_upd.R")
source("tvBAYES_log.R")
source("posterior_help.R")

# Select country of interest
# Import count data + covariates

c = "Italy"
df <- read.csv("df_bayes_2.csv") %>%
  filter(Country == c)

# Train and test (last 7 days)

policy_lag = 0
if(policy_lag != 0){
  
  l = policy_lag
  
  lag_names <- paste("lag", formatC(l, width = nchar(max(l)), flag = "0"), 
                     sep = "_")
  lag_functions <- setNames(paste("dplyr::lag(., ", l, ")"), lag_names)
  
  data_train <- df %>% 
    slice(1:(n()-7)) %>%
    mutate(Number = as.numeric(Number)) %>%
    dplyr::select(Number) %>%
    filter(row_number() > l) %>%
    pull()
  
  xreg_train = df %>%
    slice(1:(n()-7)) %>%
    dplyr::select(!c(Date,Number,Country))  %>%
    mutate_all(funs_(lag_functions)) %>%
    dplyr::select(contains("lag"))  %>%
    filter(!is.na(.)) %>%
    mutate_all(as.numeric)
  
  df_train <- df %>% 
    slice(1:(n()-7)) %>%
    filter(row_number() > l) %>%
    dplyr::select(Date, Country) %>%
    cbind(data_train, xreg_train) %>%
    rename(Number = data_train)
  
}
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
data_test <- df_test %>% 
  mutate(Number = as.numeric(Number)) %>%
  dplyr::select(Number) %>%
  pull()


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
  
  # Traceplots
  
  # traceplot(fit, "deltamufn", tv = "FALSE")
  # traceplot(fit, "deltaRfn", tv = "FALSE")
  # traceplot(fit, "pred", tv = "FALSE")
  
  # Summary of posterior distribution
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  # Posterior coefficient estimates
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=1))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At)
  
  summary <- data.frame(mut = mut$Mean[1],
                        At = At$Mean[1])
  
  summary <- data.frame(vars = colnames(summary), Coefficients=t(summary))
  rownames(summary) <- NULL
  
  # In-sample predictions
  
  vt <- predict(data=data_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 0))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  # N-ahead forecasts
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 0),
                         n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="AR", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                            model_comp = FALSE)

  
} else if (model == "ARX"){

  # Traceplots 
  
  # traceplot(fit, "deltamufn", tv = "FALSE")
  # traceplot(fit, "deltaRfn", tv = "FALSE")
  # traceplot(fit, "deltaCfn", tv = "FALSE")
  # traceplot(fit, "pred", tv = "FALSE")
  
  # Summary of posterior distribution
  
  delta <- names(fit)[names(fit) != "pred" & names(fit) != "vtp"]
  p_summary = list()
  for (j in delta){
    sum <- summary(fit, j)
    p_summary[[j]] <- sum
  }
  
  # Posterior coefficient estimates
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order_list = list(order1=1, order2=0))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order_list = list(order1 = 1, order2 = 0))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  summary <- data.frame(mut = mut$Mean[1],
                        At = At$Mean[1])
  
  Ct_mat =coeff_list[["Ct"]]
  Ct_sum = matrix(0, nrow = 1, ncol = ncol(xreg_train))
  for(k in 1:ncol(xreg_train)){
    Ct_sum[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]][1]
  }
  
  colnames(Ct_sum) <- colnames(xreg_train)
  
  summary <- cbind(summary, Ct_sum)
  summary <- data.frame(vars = colnames(summary), Coefficients=t(summary))
  rownames(summary) <- NULL
  
  # In-sample predictions
  
  vt <- predict(data=data_train, xreg = xreg_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 0))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  # N-ahead forecasts
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         xreg = xreg_train, coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 0),
                         n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="ARX", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                    model_comp = FALSE)
  

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
      
      vt <- predict(data=data_train, date = df_train %>% dplyr::select(Date),
                    coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 0))
      
      metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
      
      fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                              coeff_list=coeff_list, 
                              order_list = list(order1 = 1, order2 = 0),
                              n.ahead=7)
      
      metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
      
      fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                    model="tvAR", truth_data = df_test)
      
      metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                        model_comp = FALSE)
      
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
                   bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  At <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                    data = data_train, order_list = list(order1 = 1, order2 = 0), 
                    bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order_list = list(order1 = 1, order2 = 0),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(mut=mut, At=At, Ct=Ct)
  
  summary <- data.frame(mut=coeff_list[["mut"]][,1], At = coeff_list[["At"]][,1])
  Ct_mat =coeff_list[["Ct"]]
  Ct_sum = matrix(0, nrow = length(data_train)-1, ncol = ncol(xreg_train))
  for(k in 1:ncol(xreg_train)){
    Ct_sum[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]]
  }
  
  colnames(Ct_sum) <- colnames(xreg_train)
  summary <- cbind(summary, Ct_sum)
  
  desc_post <- describe_posterior(summary, test = c("p_direction"))
  
  vt <- predict(data=data_train, xreg = xreg_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 0))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                          xreg = xreg_train, coeff_list=coeff_list, 
                          order_list = list(order1 = 1, order2 = 0),
                          n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="tvARX", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                    model_comp = FALSE)
  
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
  
  vt <- predict(data=data_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 1))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="INGARCH", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                    model_comp = FALSE)
  
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
  
  summary <- data.frame(mut=mut[1,1], At =  At_Bt[["At"]][1,1], Bt=At_Bt[["Bt"]][1,1])
  Ct_mat =coeff_list[["Ct"]]
  Ct_sum = matrix(0, nrow = 1, ncol = ncol(xreg_train))
  for(k in 1:ncol(xreg_train)){
    Ct_sum[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]][1]
  }
  
  colnames(Ct_sum) <- colnames(xreg_train)
  summary <- cbind(summary, Ct_sum)
  summary <- data.frame(vars = colnames(summary), Coefficients=t(summary))
  rownames(summary) <- NULL
  
  # Predictions 
  
  vt <- predict(data=data_train, xreg=xreg_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 1))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  fcst <- forecast_ahead(data=data_train, xreg =xreg_train, 
                         date = df_train %>% pull(Date), coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="INGARCHX", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                    model_comp = FALSE)
  
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
  
  vt <- predict(data=data_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, 
                order_list = list(order1 = 1, order2 = 1))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  fcst <- forecast_ahead(data=data_train, date = df_train %>% pull(Date), 
                         coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="tvINGARCH", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                            model_comp = FALSE)
  
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
  
  bsplines <- bsplines_fun(data_train, order_list = list(order1 = 1, order2 = 1),
                           knot = 4, norder = 4)
  
  mut <- coeff_mut(summary=p_summary[["deltamufn"]], 
                   data= data_train, order = list(order1=1, order2=1),
                   bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  At_Bt <- coeff_At_Bt(summary_R = p_summary[["deltaRfn"]], 
                       summary_K = p_summary[["deltaKfn"]],
                       data = data_train, order = list(order1 = 1, order2 = 1), 
                       bsplines_list = list(bsplines=bsplines, knot = 4, norder=4))
  
  Ct <- coeff_Ct(summary=p_summary[["deltaCfn"]],
                 data=data_train, xreg=xreg_train, order = list(order1 = 1, order2 = 1),
                 bsplines_list = list(bsplines=bsplines, knot=4, norder=4))
  
  coeff_list <- list(sigma2lat=p_summary[["sigma2latfn"]], 
                     mut=mut, At = At_Bt[["At"]], Bt=At_Bt[["Bt"]], Ct=Ct)
  
  summary <- data.frame(mut=coeff_list[["mut"]][,1], At =coeff_list[["At"]][,1],
                        Bt=coeff_list[["Bt"]][,1])
  Ct_mat =coeff_list[["Ct"]]
  Ct_sum = matrix(0, nrow = length(data_train)-1, ncol = ncol(xreg_train))
  for(k in 1:ncol(xreg_train)){
    Ct_sum[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]]
  }
  
  colnames(Ct_sum) <- colnames(xreg_train)
  summary <- cbind(summary, Ct_sum)
  
  desc_post <- describe_posterior(summary, test = c("p_direction"))
  
  # Predictions 
  
  vt <- predict(data=data_train, xreg=xreg_train, date = df_train %>% dplyr::select(Date),
                coeff_list=coeff_list, order_list = list(order1 = 1, order2 = 1))
  
  metrics_train <- pred_eval(y_pred = vt$Mean, y_true = data_train[-1])
  
  fcst <- forecast_ahead(data=data_train, xreg =xreg_train, 
                         date = df_train %>% pull(Date), coeff_list=coeff_list, 
                         order_list = list(order1 = 1, order2 = 1),
                         n.ahead=7)
  
  metrics_test <- pred_eval(y_pred = fcst$Mean, y_true = data_test)
  
  fcst_long <- fcst_long_format(point_est = "Mean", fcst=fcst, country = c, 
                                model="tvINGARCHX", truth_data = df_test)
  
  metrics_test_scoring <- eval_fcst(fcst_long, summarise_by = c("target_variable", "model"),
                                    model_comp = FALSE)
  
}
  
  
  
  
