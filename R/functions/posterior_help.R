#' Summary of posterior distribution of raw parameters
#'
#' @description Takes the output from model functions in `BAYES_log.R` or `tvBAYES_log.R` 
#' and converts it into a dataframe with mean, median and quantiles of posterior 
#' estimates
#' @param object is a list of posterior estimates as generated model functions in `BAYES_log.R` 
#' or `tvBAYES_log.R`  
#' @param delta represents the name of the parameter of interest in the previous list `object`.
#' (`deltamufn` for intercept, `deltaRfn` for autoregressive parameter, `deltaKfn` for
#' parameter wrt dependence on conditional mean, `deltaCfn` for covariates)
#' @return dataframe with mean, median and quantiles of posterior estimates

summary <- function(object, delta = "deltamufn"){
  
  # Get delta
  if (class(delta) != "character"){
    stop("Delta needs to be of character type!")
  }
  delta <- object[[delta]]
  l <- length(delta[[1]])
  chain <- matrix(unlist(delta), l)
  
  summary_posterior <- chain %>%
    as.data.frame() %>%
    mutate(Mean = apply(., 1, mean, na.rm = TRUE),
           Median = apply(.,1, median, na.rm = TRUE),
           Quantile_0.010 = apply(., 1, quantile,0.010, na.rm = TRUE),
           Quantile_0.025 = apply(., 1, quantile,0.025, na.rm = TRUE),
           Quantile_0.050 = apply(., 1, quantile,0.050, na.rm = TRUE),
           Quantile_0.100 = apply(., 1, quantile,0.100, na.rm = TRUE),
           Quantile_0.150 = apply(., 1, quantile,0.150, na.rm = TRUE),
           Quantile_0.200 = apply(., 1, quantile,0.200, na.rm = TRUE),
           Quantile_0.250 = apply(., 1, quantile,0.250, na.rm = TRUE),
           Quantile_0.300 = apply(., 1, quantile,0.300, na.rm = TRUE),
           Quantile_0.350 = apply(., 1, quantile,0.350, na.rm = TRUE),
           Quantile_0.400 = apply(., 1, quantile,0.400, na.rm = TRUE),
           Quantile_0.450 = apply(., 1, quantile,0.450, na.rm = TRUE),
           Quantile_0.500 = apply(., 1, quantile,0.500, na.rm = TRUE),
           Quantile_0.550 = apply(., 1, quantile,0.550, na.rm = TRUE),
           Quantile_0.600 = apply(., 1, quantile,0.600, na.rm = TRUE),
           Quantile_0.650 = apply(., 1, quantile,0.650, na.rm = TRUE),
           Quantile_0.700 = apply(., 1, quantile,0.700, na.rm = TRUE),
           Quantile_0.750 = apply(., 1, quantile,0.750, na.rm = TRUE),
           Quantile_0.800 = apply(., 1, quantile,0.800, na.rm = TRUE),
           Quantile_0.850 = apply(., 1, quantile,0.850, na.rm = TRUE),
           Quantile_0.900 = apply(., 1, quantile,0.900, na.rm = TRUE),
           Quantile_0.950 = apply(., 1, quantile,0.950, na.rm = TRUE),
           Quantile_0.975 = apply(., 1, quantile,0.975, na.rm = TRUE),
           Quantile_0.990 = apply(., 1, quantile,0.990, na.rm = TRUE)) %>%
    dplyr::select(Mean, Median, starts_with("Quantile")) 
  
  return(summary_posterior)
  
}

#' Summary of posterior distribution of time-costant coefficients
#'
#' @description Takes the output from model functions in `BAYES_log.R` or `tvBAYES_log.R` 
#' and returns a list of dataframes with summary statistics of posterior distribution. 
#' @param object is a list of posterior estimates as generated model functions in `BAYES_log.R` 
#' or `tvBAYES_log.R`  
#' @param delta represents the name of the parameter of interest in the previous list `object`.
#' (`deltamufn` for intercept, `deltaRfn` for autoregressive parameter, `deltaKfn` for
#' parameter wrt dependence on conditional mean, `deltaCfn` for covariates)
#' @param cov_names is the vector with covariates' names (if present)
#' @return list of dataframes with summary statistics of posterior distribution. 

summary.fit.tc <- function(object, delta = "deltamufn", cov_names = NULL){
  
  # Get delta
  if (class(delta) != "character"){
    stop("Delta needs to be of character type!")
  }
  name <- delta
  if(!is.null(cov_names)){ 
    name <- cov_names
    }
  delta <- object[[delta]]
  l <- length(delta[[1]])
  chain <- matrix(unlist(delta), l)
  
  summary_posterior_1 <- chain %>%
    as.data.frame() %>%
    mutate(Parameter = name,
           Min = apply(., 1, mean, na.rm = TRUE),
           Quantile_0.250 = apply(., 1, quantile, 0.250, na.rm = TRUE),
           Median = apply(.,1, median, na.rm = TRUE),
           Mean = apply(., 1, mean, na.rm = TRUE),
           Quantile_0.750 = apply(., 1, quantile, 0.750, na.rm = TRUE),
           Max = apply(., 1, max, na.rm = TRUE)) %>%
    dplyr::select(Parameter,Min, Quantile_0.250, Median, Mean, Quantile_0.750, Max) 
  
  chain_2 <- as.data.frame(t(matrix(unlist(delta), l)))
  if(length(name) > 1){
  colnames(chain_2) <- cov_names
  }
  desc_post <- describe_posterior(chain_2, centrality = "median", 
                                  test = c("p_direction"), ci = 0.95) #0.89
  
  return(out <- list(summary = summary_posterior_1, 
                     desc_post = desc_post))
  
}

#' Summary of posterior distribution of time-varying coefficients
#'
#' @description Takes the output from model functions in `BAYES_log.R` or `tvBAYES_log.R` 
#' and returns a list of dataframes with summary statistics of posterior distribution. 
#' @param post_coeff represents dataframe of the posterior coefficient estimate 
#' of interest as returned from `coeff_mut`, `coeff_At_Bt`, `coeff_Ct`
#' @return list of dataframes with summary statistics of posterior distribution. 
  
summary.fit.tv <- function(post_coeff){
  
  par <- post_coeff %>%
    dplyr::select(Mean)
  
  summary_posterior_1 <- par  %>%
    as.data.frame() %>%
    mutate(Min = apply(., 2, min, na.rm = TRUE),
           Quantile_0.250 = apply(., 2, quantile, 0.250, na.rm = TRUE),
           Median = apply(.,2, median, na.rm = TRUE),
           Mean = apply(., 2, mean, na.rm = TRUE),
           Quantile_0.750 = apply(., 2, quantile, 0.750, na.rm = TRUE),
           Max = apply(., 2, max, na.rm = TRUE)) %>%
    dplyr::select(Min, Quantile_0.250, Median, Mean, Quantile_0.750, Max) %>%
    head(1)
  
  desc_post <- describe_posterior(par, centrality = "median", 
                                  test = c("p_direction"), ci = 0.95)
  
  return(out <- list(summary = summary_posterior_1, 
                     desc_post = desc_post))
  
}

#' B-splines
#'
#' @description Takes input data and details of model specifications and returns
#' a dataframe of b-splines to compute coefficient estimates
#' @param data  represents count data used in model fit
#' @param order_list is a list of the specification of model order (p,q)
#' @param knot is the number of equidistant knots for B-spline
#' @param norder is the order of B-splines
#' @return dataframe of b-splines

bsplines_fun <- function(data, order_list = list(order1 = 1, order2 = 1),
                     knot = 4, norder = 4){
  
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
  } else{
    order = order_list[["order1"]]
  } 
  
  #B-splines
  time <- (1:length(data)) / length(data)
  J       <- knot + norder-1
  timesp  <- bsplineS(time,  breaks=seq(0,1,1/knot), norder = norder)
  timespI <- timesp
  if(order>0){timespI <- timesp[-(1:order), ]}
  
  return(timespI)
}
  
#' Posterior coefficient estimates
#'
#' @description Takes the summary of posterior distribution, data, details of 
#' model specification (order p, q and b-splines) and returns 
#' a dataframe of posterior coefficient estimates
#' @param summary is the output from `summary` function
#' @param data represents count data used in model fit
#' @param order_list is a list of the specification of model order (p,q)
#' @param bsplines_list is a list of the specification of b-splines
#' @return dataframe of posterior coefficient estimates

coeff_mut <- function(summary, data, order_list = list(order1 = 1, order2 = 1), 
                                   bsplines_list = NULL){
  
 
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
  } else{
    order = order_list[["order1"]]
  } 
  
  mut = matrix(0,length(data)-order, ncol(summary))
  
  if(!is.null(bsplines_list)){
    for (i in 1:ncol(summary)){
      mut[,i] <- bsplines_list[["bsplines"]] %*% summary[,i]
    }
    colnames(mut) <- colnames(summary)
  } else {
    for (i in 1:ncol(summary)){
      mut[,i] <- rep(summary[,i], (length(data)-order))
    }
    colnames(mut) <- colnames(summary)
  }
  
  return(data.frame(mut))
}


coeff_At_Bt <- function(summary_R, summary_K = NULL, data, 
                           order_list = list(order1 = 1, order2 = 0), 
                           bsplines_list = NULL){
  
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
  
    At_mat <- matrix(0, length(data) - order, ncol(summary_R))
    Bt_mat <- matrix(0, length(data) - order, ncol(summary_K))
    
    if(!is.null(bsplines_list)){
      J = bsplines_list[["knot"]] + bsplines_list[["norder"]]-1
      for (i in 1:ncol(summary_R)){
        
        At <- matrix(0, length(data) - order, order1)
        # deltaA <- (summary_R[,i]+summary_K[,i])/2
        deltaA <- (summary_R[,i])
        for(j in 1:order1){
          At[, j] <- bsplines_list[["bsplines"]] %*% deltaA[(j-1)*J+1:J]
        }
        At_mat[,i] <- At
        
        Bt <- matrix(0, length(data) - order, order2)
        # deltaB <- (summary_R[,i]-summary_K[,i])/2
        deltaB <- (summary_K[,i])
        for(j in 1:order2){
          Bt[, j] <- bsplines_list[["bsplines"]] %*% deltaB[(j-1)*J+1:J]
        }
        Bt_mat[,i] <- Bt
        
      }
      colnames(At_mat) <- colnames(summary_R)
      colnames(Bt_mat) <- colnames(summary_K)
    } else {
      for (i in 1:ncol(summary_R)){
        
        At <- matrix(0, length(data) - order, order1)
        # deltaA <- (summary_R[,i]+summary_K[,i])/2
        deltaA <- (summary_R[,i])
        for(j in 1:order1){
          At[, j] <- deltaA[j]
        }
        At_mat[,i] <- At
        
        Bt <- matrix(0, length(data) - order, order2)
        # deltaB <- (summary_R[,i]-summary_K[,i])/2
        deltaB <- (summary_K[,i])
        for(j in 1:order2){
          Bt[, j] <- deltaB[j]
        }
        Bt_mat[,i] <- Bt
        
      }
      colnames(At_mat) <- colnames(summary_R)
      colnames(Bt_mat) <- colnames(summary_K)
    }
    } else{
      
      order = order_list[["order1"]]
    At_mat <- matrix(0, length(data) - order, ncol(summary_R))
    
    if (!is.null(bsplines_list)){
      
      J = bsplines_list[["knot"]] + bsplines_list[["norder"]]-1
    for (i in 1:ncol(summary_R)){
      
      At <- matrix(0, length(data) - order, order)
      deltaA <- (summary_R[,i])
      for(j in 1:order){
        At[,j] <- bsplines_list[["bsplines"]] %*% deltaA[(j-1)*J+1:J]
      }
      At_mat[,i] <- At
      
    }
    colnames(At_mat) <- colnames(summary_R)
    
    } else {
      
      for (i in 1:ncol(summary_R)){
        
        At <- matrix(0, length(data) - order, order)
        deltaA <- (summary_R[,i])
        for(j in 1:order){
          At[, j] <- deltaA[j]
        }
        At_mat[,i] <- At
        
      }
      colnames(At_mat) <- colnames(summary_R)
    }
    }
  if(order_list[["order2"]] != 0){ out <- list(At = data.frame(At_mat), Bt = data.frame(Bt_mat))}
  if(order_list[["order2"]] == 0){ out <- data.frame(At_mat)}
  return(out)
}
  

coeff_Ct <- function(summary, data, xreg, 
                       order_list = list(order1 = 1, order2 = 1), 
                       bsplines_list = NULL){
  
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
  } else{
    order = order_list[["order1"]]
  } 
  
  if(!is.null(bsplines_list)){
    J = bsplines_list[["knot"]] + bsplines_list[["norder"]]-1
    n = ncol(xreg)
    Ct_com <- list()
    Ct_mat <- matrix(0, length(data) - order, ncol(summary))
    for(k in 1:n){
      for (i in 1:ncol(summary)){
        Ct <- matrix(0, length(data) - order, 1)
        deltaC <- summary[,i]
        Ct[, 1] <- bsplines_list[["bsplines"]]  %*% deltaC[(k-1)*J+1:J]
        Ct_mat[,i] <- Ct
      }
      colnames(Ct_mat) <- colnames(summary)
      Ct_com[[paste0("Ct", k)]] <- data.frame(Ct_mat)
  }} else{
    n = ncol(xreg)
    Ct_com <- list()
    Ct_mat <- matrix(0, length(data) - order, ncol(summary))
    
    for(k in 1:n){
      for (i in 1:ncol(summary)){
        Ct <- matrix(0, length(data) - order, 1)
        deltaC <- summary[,i]
        Ct[, 1] <- deltaC[k]
        Ct_mat[,i] <- Ct
      }
      colnames(Ct_mat) <- colnames(summary)
      Ct_com[[paste0("Ct", k)]] <- data.frame(Ct_mat)
    }
    }
  return(Ct_com)
}


#' In-sample predictions 
#'
#' @description Takes data, coefficient list and specification of model order and 
#' returns a dataframe of the distribution of in-sample predictions
#' @param data represents count data used in model fit
#' @param date represents the date object corresponding to count data used in model fit
#' @param xreg is the dataframe containing values for covarates used in model fit
#' @param coeff_list is the list of posterior coefficients
#' @param order_list is a list of the specification of model order (p,q)
#' @return dataframe of the distribution of in-sample predictions

predict <- function(data, date, xreg = NULL, coeff_list,
                       order_list = list(order1 = 1, order2 = 0)){
  
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
    
    Y <- data 
    if(order>0){Y <- data[-(1:order)]}
    X <- NULL
    if(order > 0){
      for(j in 1:order){
        ind <- (order-j+1):(length(data) - j)
        X <- cbind(X, array(data[ind]))
      } 
    }
    date <- date %>% filter(row_number() > order)
    
    temp <- as.data.frame(coeff_list[["sigma2lat"]])
    mut <- as.data.frame(coeff_list[["mut"]])
    At <- as.data.frame(coeff_list[["At"]])
    Bt <- as.data.frame(coeff_list[["Bt"]])
    
    if(is.null(xreg)){
      
      vt_mat = matrix(0, nrow=length(data)-order, ncol=ncol(mut))
      for (i in 1:ncol(mut)){
        temp_s <- temp[,i]
        mut_s <- mut[,i]
        At_s <- At[,i]
        Bt_s <- Bt[,i]
        vt = matrix(0, nrow=length(data)-order)
        for (j in 1:(length(data)-order)){
          vt[j,] <- mut_s[j] + sum(At_s[j]*log1p(X[j,]))+ sum(Bt_s[j]*temp_s[(order2+j-1):j])
          temp_s    <- c(temp_s, vt[j,])
        }
        vt_mat[,i] <- vt
      colnames(vt_mat) <- colnames(mut)
      
      }
      } else{
      
      Ct_mat <- coeff_list[["Ct"]]
      n <- ncol(xreg)
      if(order>0){xreg <- xreg[-(1:order),]}
      vt_mat = matrix(0, nrow=length(data)-order, ncol=ncol(mut))
      Ct_s = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
      for (i in 1:ncol(mut)){
        temp_s <- temp[,i]
        mut_s <- mut[,i]
        At_s <- At[,i]
        Bt_s <- Bt[,i]
      for(k in 1:n){
        Ct_s[,k] <- Ct_mat[[paste0("Ct", k)]][,i]
      }
        vt = matrix(0, nrow=length(data)-order)
        for (j in 1:(length(data)-order)){
          vt[j,] <- mut_s[j] + sum(At_s[j]*log1p(X[j,]))+ sum(Bt_s[j]*temp_s[(order2+j-1):j])+ sum(Ct_s[j,]*xreg[j,])
          temp_s    <- c(temp_s, vt[j,])
        }
        vt_mat[,i] <- vt
        colnames(vt_mat) <- colnames(mut)
      }
      }
    } else{
        
        order = order_list[["order1"]]
        Y <- data 
        if(order>0){Y <- data[-(1:order)]}
        X <- NULL
        if(order > 0){
          for(j in 1:order){
            ind <- (order-j+1):(length(data) - j)
            X <- cbind(X, array(data[ind]))
          } 
        }
        
        date <- date %>% filter(row_number() > order)
    
    mut <- as.data.frame(coeff_list[["mut"]])
    At <- as.data.frame(coeff_list[["At"]])

  if(is.null(xreg)){
    
    vt_mat = matrix(0, nrow=length(data)-order, ncol=ncol(mut))
    for (i in 1:ncol(mut)){
      comp2 = array(At[,i] * log1p(X))
      vt_mat[,i] <- array(mut[,i]) + array(comp2) 
    }
    colnames(vt_mat) <- colnames(mut)
    
  } else{
    
    Ct_mat <- coeff_list[["Ct"]]
    n <- ncol(xreg)
    if(order>0){xreg <- xreg[-(1:order),]}
    vt_mat = matrix(0, nrow=length(data)-order, ncol=ncol(mut))
    Ct = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
      for (i in 1:ncol(mut)){
        for(k in 1:n){
          Ct[,k] <- Ct_mat[[paste0("Ct", k)]][,i]
        }
        comp2 = array(At[,i] * log1p(X))
        comp3 <- array(rowSums(Ct * xreg))
        vt_mat[,i] <- array(mut[,i]) + array(comp2) + array(comp3)
      }
      colnames(vt_mat) <- colnames(mut)
    }
    }
  
  vt_mat <- data.frame(cbind(date, exp(vt_mat)))
  return(vt_mat)
}

#' Evaluation metrics
#'
#' @description Takes true and predicted data and return a dataframe with MLmetrics
#' wrt model prediction error. 
#' @param y_pred represents predicted data
#' @param y_true represents count data used in model fit
#' @return dataframe with MLmetrics wrt model prediction error. 

pred_eval <- function(y_pred, y_true){
  
  MAE = MLmetrics::MAE(y_pred, y_true)
  MAPE = MLmetrics::MAPE(y_pred, y_true)
  MedianAE = MLmetrics::MedianAE(y_pred, y_true)
  MedianAPE = MLmetrics::MedianAPE(y_pred, y_true)
  MSE = MLmetrics::MSE(y_pred, y_true)
  RMSE = MLmetrics::RMSE(y_pred, y_true)
  
  pred_metrics <- data.frame(MAE,MAPE,MedianAE,MedianAPE,MSE,RMSE)
  
  return(pred_metrics)
}

#' N-step ahead forecasts
#'
#' @description Takes data, coefficient list and specification of model orders and 
#' return n-step ahead forecasts
#' @param data represents count data used in model fit
#' @param xreg_train is the dataframe of covariates used in model fit
#' @param newxreg represents the dataframe containing new values for the covariates 
#' to be used for prediction. If newxreg is omitted the last known values of the covariates 
#' are used for prediction
#' @param date represents the date object corresponding to count data used in model fit
#' @param coeff_list is the list of posterior coefficients
#' @param order_list is a list of the specification of model order (p,q)
#' @param validation logical value specifying whether predictions will be validated
#' on a test set (validation = TRUE) or not (validation = FALSE)
#' @param n.ahead positive integer value giving the number of steps ahead for which 
#' predictions should be made
#' @return dataframe of distribution of n-step ahead forecasts

forecast_ahead <- function(data, xreg_train = NULL, newxreg = NULL, date, coeff_list = list(),
                           order_list = list(order1 = 1, order2 = 1),
                           validation = TRUE, n.ahead = 1){
  
  lastDate <- date %>% last() %>% as.Date()
  if (order_list[["order2"]] != 0){
    order1 = order_list[["order1"]]
    order2 = order_list[["order2"]]
    order = max(order1, order2)
    
    Y <- data
    if(order>0){Y <- data[-(1:order)]}
    X <- NULL
    if(order > 0){
      for(j in 1:order){
        ind <- (order-j+1):(length(data) - j)
        X <- cbind(X, array(data[ind]))
      } 
    }
    
    temp <- as.data.frame(coeff_list[["sigma2lat"]])
    mut <- as.data.frame(coeff_list[["mut"]])
    At <- as.data.frame(coeff_list[["At"]])
    Bt <- as.data.frame(coeff_list[["Bt"]])
    
    if(is.null(xreg_train)){
      
      vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
      pred <- data.frame(vt); colnames(pred)=colnames(mut)
      
      for(i in 1:ncol(mut)){
        temp_s <- temp[,i]
        mut_s <- data.frame(mut[,i])
        At_s <- data.frame(At[,i])
        Bt_s <- data.frame(Bt[,i])
        data_pred <- data 
        for (k in 1:nrow(pred)){
          
          temp_s <- temp[,i]
          X <- NULL
          if(order > 0){
            for(j in order){
              ind <- (order-j+1):(length(data_pred) - j)
              X <- cbind(X, array(data_pred[ind]))
            } 
          }
          
          mut_k = mut_s %>%
            slice(rep(n(), each = 1))
          mut_s = rbind(mut_s, mut_k)
          
          At_k = At_s %>% 
            slice(rep(n(), each = 1))
          At_s = rbind(At_s,At_k)
          
          Bt_k = Bt_s %>% 
            slice(rep(n(), each = 1))
          Bt_s = rbind(Bt_s,Bt_k)
          
          vt <- matrix(rep(0, (length(data_pred)-order)))
          for(j in 1:(length(data_pred)-order)){
            vt[j,] <- mut_s[j,] + sum(At_s[j,]*log1p(X[j, ]))+ sum(Bt_s[j,]*temp_s[(order2+j-1):j])
            temp_s    <- c(temp_s, vt[j,])
          } 
          
          t = length(temp_s)
          pred[k,i] <- as.numeric(exp(mut_s[t,] + sum(At_s[t,]*log1p(data_pred[t]))+ sum(Bt_s[t,]*temp_s[(t)])))
          
          data_pred <- c(data_pred, pred[k,i])
          
        }
      }
      p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
      pred <- cbind(p, pred)
    } else{
      if (validation == FALSE){
        Ct_mat <- coeff_list[["Ct"]]
        n <- ncol(xreg_train)
        if(order>0){xreg <- xreg_train[-(1:order),]}
        vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
        pred <- data.frame(vt); colnames(pred)=colnames(mut)
        # Ct_s = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
        
        for(i in 1:ncol(mut)){
          temp_s <- temp[,i]
          mut_s <- data.frame(mut[,i])
          At_s <- data.frame(At[,i])
          Bt_s <- data.frame(Bt[,i])
          data_pred <- data 
          xreg_pred <- xreg
          Ct_s = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
          for(s in 1:n){
            Ct_s[,s] <- Ct_mat[[paste0("Ct", s)]][,i]
          }
          Ct_s <- data.frame(Ct_s)
          for (k in 1:nrow(pred)){
            
            temp_s <- temp[,i]
            xreg_pred_k = xreg_pred %>%
              slice(rep(n(), each = (k)))
            xreg_pred = rbind(xreg_pred, xreg_pred_k)
            
            X <- NULL
            if(order > 0){
              for(j in order){
                ind <- (order-j+1):(length(data_pred) - j)
                X <- cbind(X, array(data_pred[ind]))
              } 
            }
            
            mut_k = mut_s %>%
              slice(rep(n(), each = 1))
            mut_s = rbind(mut_s, mut_k)
            
            At_k = At_s %>% 
              slice(rep(n(), each = 1))
            At_s = rbind(At_s,At_k)
            
            Bt_k = Bt_s %>% 
              slice(rep(n(), each = 1))
            Bt_s = rbind(Bt_s,Bt_k)
            
            Ct_k = Ct_s %>% 
              slice(rep(n(), each = 1))
            Ct_s = rbind(Ct_s,Ct_k)
            
            vt <- matrix(rep(0, (length(data_pred)-order)))
            for(j in 1:(length(data_pred)-order)){
              vt[j,] <- mut_s[j,] + sum(At_s[j, ]*log1p(X[j, ]))+ sum(Bt_s[j, ]*temp_s[(order2+j-1):j]) + sum(Ct_s[j,]*xreg_pred[j,])
              temp_s    <- c(temp_s, vt[j,])
            } 
            
            t = length(temp_s)
            pred[k,i] <- as.numeric(exp(mut_s[t,] + sum(At_s[t, ]*log1p(data_pred[t]))+ sum(Bt_s[t, ]*temp_s[(t)]) + sum(Ct_s[t,]*xreg_pred[t,])))
            
            data_pred <- c(data_pred, pred[k,i])
            
          }
        }
        p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
        pred <- cbind(p, pred)
      } else if(validation == TRUE){
        
        Ct_mat <- coeff_list[["Ct"]]
        n <- ncol(xreg_train)
        if(order>0){xreg <- xreg_train[-(1:order),]}
        vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
        pred <- data.frame(vt); colnames(pred)=colnames(mut)
        # Ct_s = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
        
        for(i in 1:ncol(mut)){
          temp_s <- temp[,i]
          mut_s <- data.frame(mut[,i])
          At_s <- data.frame(At[,i])
          Bt_s <- data.frame(Bt[,i])
          data_pred <- data 
          xreg_pred <- xreg
          Ct_s = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
          for(s in 1:n){
            Ct_s[,s] <- Ct_mat[[paste0("Ct", s)]][,i]
          }
          Ct_s <- data.frame(Ct_s)
          for (k in 1:nrow(pred)){
            
            temp_s <- temp[,i]
            xreg_pred_k = newxreg %>%
              filter(row_number()==k)
            xreg_pred = rbind(xreg_pred, xreg_pred_k)
            
            X <- NULL
            if(order > 0){
              for(j in order){
                ind <- (order-j+1):(length(data_pred) - j)
                X <- cbind(X, array(data_pred[ind]))
              } 
            }
            
            mut_k = mut_s %>%
              slice(rep(n(), each = 1))
            mut_s = rbind(mut_s, mut_k)
            
            At_k = At_s %>% 
              slice(rep(n(), each = 1))
            At_s = rbind(At_s,At_k)
            
            Bt_k = Bt_s %>% 
              slice(rep(n(), each = 1))
            Bt_s = rbind(Bt_s,Bt_k)
            
            Ct_k = Ct_s %>% 
              slice(rep(n(), each = 1))
            Ct_s = rbind(Ct_s,Ct_k)
            
            vt <- matrix(rep(0, (length(data_pred)-order)))
            for(j in 1:(length(data_pred)-order)){
              vt[j,] <- mut_s[j,] + sum(At_s[j, ]*log1p(X[j, ]))+ sum(Bt_s[j, ]*temp_s[(order2+j-1):j]) + sum(Ct_s[j,]*xreg_pred[j,])
              temp_s    <- c(temp_s, vt[j,])
            } 
            
            t = length(temp_s)
            pred[k,i] <- as.numeric(exp(mut_s[t,] + sum(At_s[t, ]*log1p(data_pred[t]))+ sum(Bt_s[t, ]*temp_s[(t)]) + sum(Ct_s[t,]*xreg_pred[t,])))
            
            data_pred <- c(data_pred, pred[k,i])
            
          }
        }
        p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
        pred <- cbind(p, pred)
      }
    }
  }else{
    
    order = order_list[["order1"]]
    Y <- data 
    if(order>0){Y <- data[-(1:order)]}
    X <- NULL
    if(order > 0){
      for(j in 1:order){
        ind <- (order-j+1):(length(data) - j)
        X <- cbind(X, array(data[ind]))
      } 
    }
    
    mut <- data.frame(coeff_list[["mut"]])
    At <- data.frame(coeff_list[["At"]])
    
    if(is.null(xreg_train)){
      
      vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
      pred <- data.frame(vt); colnames(pred)=colnames(mut)
      
      for(i in 1:ncol(mut)){
        mut_k = data.frame(mut[,i]) %>% filter(row_number()==n()) %>% pull()
        At_k = data.frame(At[,i]) %>% filter(row_number()==n()) %>% pull()
        data_pred <- data.frame(data) %>% filter(row_number()==n()) %>% pull()
        for (k in 1:nrow(pred)){
          
          # t = length(data_pred)
          pred[k,i] <- as.numeric(exp(array(mut_k) + array(At_k* log1p(data_pred[k]))))
          data_pred <- c(data_pred, pred[k,i])
          
        }
      }
      
      p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
      pred <- cbind(p, pred)
    } else{
      if (validation == FALSE){
        Ct_mat <- coeff_list[["Ct"]]
        n <- ncol(xreg_train)
        if(order>0){xreg <- xreg_train[-(1:order),]}
        vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
        pred <- data.frame(vt); colnames(pred)=colnames(mut)
        Ct = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
        
        for(i in 1:ncol(mut)){
          mut_k = data.frame(mut[,i]) %>% filter(row_number()==n()) %>% pull()
          At_k = data.frame(At[,i]) %>% filter(row_number()==n()) %>% pull()
          for(s in 1:n){
            Ct[,s] <- Ct_mat[[paste0("Ct", s)]][,i]
          }
          Ct_k = data.frame(Ct) %>% filter(row_number()==n()) 
          data_pred <- data.frame(data) %>% filter(row_number()==n()) %>% pull()
          xreg_pred <- data.frame(xreg) %>% filter(row_number()==n()) 
          for (k in 1:nrow(pred)){
            
            pred[k,i] <- as.numeric(exp(array(mut_k) + array(At_k* log1p(data_pred[k])) + array(rowSums(Ct_k * xreg_pred))))
            data_pred <- c(data_pred, pred[k,i])
            
          }
        }
        p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
        pred <- cbind(p, pred)
      } else if (validation == TRUE){
        
        Ct_mat <- coeff_list[["Ct"]]
        n <- ncol(newxreg)
        # if(order>0){xreg <- xreg_train[-(1:order),]}
        vt <- matrix(0, nrow=n.ahead, ncol=ncol(mut))
        pred <- data.frame(vt); colnames(pred)=colnames(mut)
        Ct = matrix(0, nrow = length(data)-order, ncol = ncol(newxreg))
        
        for(i in 1:ncol(mut)){
          mut_k = data.frame(mut[,i]) %>% filter(row_number()==n()) %>% pull()
          At_k = data.frame(At[,i]) %>% filter(row_number()==n()) %>% pull()
          for(s in 1:n){
            Ct[,s] <- Ct_mat[[paste0("Ct", s)]][,i]
          }
          Ct_k = data.frame(Ct) %>% filter(row_number()==n()) 
          data_pred <- data.frame(data) %>% filter(row_number()==n()) %>% pull()
          # xreg_pred <- data.frame(xreg) %>% filter(row_number()==n()) 
          for (k in 1:nrow(pred)){
            
            xreg_pred <- data.frame(newxreg) %>% filter(row_number()==k)
            pred[k,i] <- as.numeric(exp(array(mut_k) + array(At_k* log1p(data_pred[k])) + array(rowSums(Ct_k * xreg_pred))))
            data_pred <- c(data_pred, pred[k,i])
            
          }
        }
        p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
        pred <- cbind(p, pred)
        
      }
    }
  }
  return(pred)
}

extract_numeric <- function(x) {
  as.numeric(gsub("[^0-9.-]+", "", as.character(x)))
}

#' Long format of model predictions
#'
#' @description Takes model predictions and returns a long dataframe with 
#' point estimate and quantiles of predicted values and observed ones 
#' @param point_est represents the point estimate to summarise predictive distribution
#' @param country is the character value with the name of the country of interest
#' @param model represents the character value specifying model of interest
#' @param truth_data is the dataframe of true data
#' @return a long dataframe with point estimate and quantiles of predicted values 
#' and observed ones 

fcst_long_format <- function(point_est = "Mean", fcst, country, model, truth_data){
  
  fcst_tot = c()
  for (i in 1:nrow(fcst)){
    
    fcst_2 <- fcst %>%
      slice(i) %>%
      dplyr::select(all_of(point_est), starts_with("Quantile")) %>%
      gather() %>%
      mutate(quantile = extract_numeric(key),
             key = ifelse(key == "Mean", "point", key),
             key = ifelse(str_detect(key, "Quantile"), "quantile", key),
             horizon = i,
             temporal_resolution = "day", 
             target_end_date = fcst$forecast_date[i]) %>%
      rename(type = key) %>%
      dplyr::select(quantile, value, type, horizon, temporal_resolution, target_end_date) 
    
    fcst_tot <- fcst_tot %>%
      rbind(fcst_2)
  }
  
  fcst_tot <- fcst_tot %>%
    mutate(forecast_date = fcst$forecast_date[1]-1,
           location = country, 
           target_variable = "inc cases", 
           model = model) %>%
    rename(prediction = value) %>%
    dplyr::select(model, forecast_date, target_variable, 
                  target_end_date, location, quantile, 
                  prediction, type)
  
  truth <- truth_data %>%
    mutate(target_end_date= as.Date(Date)) %>%
    rename(true_value = Number) %>%
    dplyr::select(target_end_date, true_value)
  
  data <- scoringutils::merge_pred_and_obs(fcst_tot, truth,
                                           join = "full")
  
  return(data)
}

#' Scoringutils evaluation metrics
#'
#' @description Takes the long format of model predictions and returns a collection 
#' of different metrics and scoring rules as in `scoringutils` package.
#' @param fcst_long_format represents a data.frame with the predictions and observations
#' as returned from `fcst_long_format` function
#' @param summarise_by is the character vector of columns to group the summary by. 
#' @param model_comp logical, whether or not to compute relative performance between models. 
#' @return dataframe with a collection of different metrics and scoring rules as 
#' in `scoringutils` package.

eval_fcst <- function(fcst_long_format, summarise_by = c("target_variable", "model"),
                      model_comp = FALSE){
  
  if (model_comp == TRUE){
  
    coverage <- eval_forecasts(
      fcst_long_format,
    summarise_by = c(summarise_by, "range"),
    compute_relative_skill = TRUE) %>%
    dplyr::filter(range %in% c(50, 95)) %>%
    dplyr::select(model, coverage, range) %>%
    tidyr::pivot_wider(names_from = range, values_from = coverage,
                       names_prefix = "Coverage ")

  table <- eval_forecasts(fcst_long_format, summarise_by = summarise_by,
                          compute_relative_skill = TRUE) %>%
    dplyr::left_join(coverage, by = "model")
  
  # table <- eval_forecasts(fcst_long_format, summarise_by = summarise_by, 
  #                         compute_relative_skill = TRUE) 
  
  } else{
    
    coverage <- eval_forecasts(
      fcst_long_format,
      summarise_by = c(summarise_by, "range"),
      compute_relative_skill = FALSE) %>%
      dplyr::filter(range %in% c(50, 95)) %>%
      dplyr::select(model, coverage, range) %>%
      tidyr::pivot_wider(names_from = range, values_from = coverage,
                         names_prefix = "Coverage ")

    table <- eval_forecasts(fcst_long_format, summarise_by = summarise_by,
                            compute_relative_skill = FALSE) %>%
      dplyr::left_join(coverage, by = "model")
    
    # table <- eval_forecasts(fcst_long_format, summarise_by = summarise_by, 
    #                         compute_relative_skill = FALSE) 
    
  }
  
}

#' COVID-19 European forecast hub (weekly forecasts)
#'
#' @description Takes the wide format of model predictions and returns a dataframe
#' of long predictions aggregated by epidemiological week (as required by COVID-19 
#' European Forecast Hub)
#' @param point_est represents the point estimate to summarise predictive distribution
#' @param fcst represents dataframe of model predictions
#' @param country is the character value with the name of the country of interest
#' @return a dataframe of long predictions aggregated by epidemiological week 

eu_fcst_hub <- function(point_est = "Mean", fcst, country){
  
  fcst_tot = c()
  for (i in 1:nrow(fcst)){
    
    fcst_2 <- fcst %>%
      slice(i) %>%
      dplyr::select(all_of(point_est), starts_with("Quantile")) %>%
      gather() %>%
      mutate(quantile = extract_numeric(key),
             key = ifelse(key == "Mean", "point", key),
             key = ifelse(str_detect(key, "Quantile"), "quantile", key),
             horizon = i,
             temporal_resolution = "day", 
             target_end_date = fcst$forecast_date[i]) %>%
      rename(type = key) %>%
      dplyr::select(quantile, value, type, horizon, temporal_resolution, target_end_date) 
    
    fcst_tot <- fcst_tot %>%
      rbind(fcst_2)
  }
  
  fcst_tot <- fcst_tot %>%
    mutate(epi_week = epiweek(target_end_date), 
           week_day = weekdays(target_end_date))
  
  hub <- fcst_tot %>%
    group_by(quantile, epi_week) %>%
    summarise(value = sum(value)) %>%
    ungroup()
  
  week <- fcst_tot %>%
    filter(week_day == "Saturday") %>%
    dplyr::select(target_end_date, epi_week) %>%
    distinct(epi_week, target_end_date)
  
  hub <- plyr::join(hub, week, by = "epi_week")
  
  hub <- hub %>%
    mutate(type = ifelse(is.na(quantile), "point", "quantile"))
  
  report_date <- fcst %>%
    dplyr::select(forecast_date) %>%
    mutate(day = weekdays(forecast_date)) %>%
    filter(day == "Monday") %>%
    filter(row_number() == 1) %>%
    pull(forecast_date)
  
  library(countrycode)
  location_code = countrycode(country, origin = 'country.name', destination = 'iso2c')
  
  hub <- hub %>%
    mutate(forecast_date = report_date,
           week_ahd = seplyr::add_group_indices(., c("target_end_date"), "group_ID")$group_ID,
           target = paste(week_ahd, "wk ahead inc case"),
           location = location_code, 
           scenario_id = "forecast") %>%
    dplyr::select(forecast_date, target, target_end_date, location, quantile, value,
                  type, scenario_id) 
  
  return(hub)
  
}
  
