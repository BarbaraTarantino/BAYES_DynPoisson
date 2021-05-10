# Traceplot

traceplot <- function(object, delta, thin = 0, tv = "FALSE"){
  
  # Get delta
  if (class(delta) != "character"){
    stop("Delta needs to be of character type!")
  }
  par <- object[[delta]]
  chain  <- matrix(unlist(par), length(par[[1]]))
  
  # Thin
  if (thin > 0){
    colnum <- ncol(chain)
    thin <- abs(round(thin))
    if(thin > colnum) stop("Thin exceeds number of MCMC samples.")
    keepcols <- which(rep(1:thin, len=colnum) == thin)
    chain <- t(chain[,keepcols])
  }
  
  # Traceplot
  if (tv == "FALSE"){
    if(delta == "sigma2latfn" | delta == "pred" ){
      n = ncol(chain)
      return(plot(chain[2:n] - chain[1:(n-1)])^2)
    } else{
      n = ncol(chain)
      return(plot(chain[, 2:n] - chain[, 1:(n-1)])^2)
      }} else if (tv == "TRUE"){
    if(delta == "sigma2latfn" | delta == "pred"){
      n = ncol(chain)
      return(plot(chain[2:n] - chain[1:(n-1)])^2)
    } else{
      n = ncol(chain)
      return(plot(colMeans((chain[, 2:n] - chain[, 1:(n-1)])^2)))
  }}
  
}


# Summary of posterior distribution (mean, median and quantiles)

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
  
# B-splines

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
  
# Posterior coefficient estimates (mut, At, Bt, Ct)

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


# In-sample Predictions

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

# N-step ahead forecast

forecast_ahead <- function(data, xreg = NULL, date, coeff_list = list(),
                              order_list = list(order1 = 1, order2 = 1), n.ahead = 1){
  
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
    
    if(is.null(xreg)){
      
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
      
      Ct_mat <- coeff_list[["Ct"]]
      n <- ncol(xreg)
      if(order>0){xreg <- xreg[-(1:order),]}
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
    
    mut <- data.frame(coeff_list[["mut"]])
    At <- data.frame(coeff_list[["At"]])
    
    if(is.null(xreg)){
      
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
      
      Ct_mat <- coeff_list[["Ct"]]
      n <- ncol(xreg)
      if(order>0){xreg <- xreg[-(1:order),]}
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
    }
  }
  return(pred)
}

fcst_plot <- function(point_est = "Mean", fcst, country, model, truth_data){
  
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
  
  truth_data <- truth_data %>%
    dplyr::select(Country, Date,Number) %>%
    mutate(Date = as.Date(Date),
           Number = as.numeric(Number)) %>%
    rename(target_end_date = Date,
           location = Country,
           true_value = Number) %>%
    mutate(target_variable = unique(fcst_tot$target_variable)) %>%
    dplyr::select(target_end_date, true_value, 
                  target_variable, location)
  
  data_plot <- truth_data %>%
    dplyr::bind_rows(fcst_tot) %>%
    filter(target_end_date > "2021-02-28")
  
  plot <- scoringutils::plot_predictions(data_plot,
                                         x = "target_end_date",
                                         facet_formula = ~ model,
                                         ncol = 3,
                                         # facet_formula = model ~ target_variable + loc,
                                         # facet_wrap_or_grid = "facet",
                                         allow_truth_without_pred = FALSE,
                                         scales = "free") + 
    # ggplot2::ggtitle(paste0("Predictions for incident ", target_variable,  "s")) + 
    ggplot2::theme(legend.position = "bottom", 
                   strip.placement = "outside") + 
    scale_y_continuous(labels = scales::comma) + 
    expand_limits(y = 0) + 
    coord_cartesian(ylim = c(0, NA))
  
  return(print(plot))
}

extract_numeric <- function(x) {
  as.numeric(gsub("[^0-9.-]+", "", as.character(x)))
}

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
  
  hub <- hub %>%
    mutate(forecast_date = report_date,
           week_ahd = seplyr::add_group_indices(., c("target_end_date"), "group_ID")$group_ID,
           target = paste(week_ahd, "wk ahead inc case"),
           location = "IT", 
           scenario_id = "forecast") %>%
    dplyr::select(forecast_date, target, target_end_date, location, quantile, value,
                  type, scenario_id) 
  
  return(hub)
  
}
  

post_contribution_plot <- function(data, xreg = NULL, coeff_list,
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
    
    xreg <- xreg[-(1:order),]
    n <- ncol(xreg)
    
    temp <- as.numeric(as.data.frame(coeff_list[["sigma2lat"]][["Mean"]]))
    mut <- as.data.frame(coeff_list[["mut"]][["Mean"]])
    At <- as.data.frame(coeff_list[["At"]][["Mean"]])
    Bt <- as.data.frame(coeff_list[["Bt"]][["Mean"]])
    Ct_mat <- coeff_list[["Ct"]]
    Ct = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
    for(k in 1:n){
      Ct[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]]
    }
    Ct <- as.data.frame(Ct)
    
    Int <-c(); AR <- c()
    MA <- c(); XREG = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
    
    vt = matrix(0, nrow=length(data)-order)
    for (j in 1:(length(data)-order)){
      vt[j] <- mut[j,] + sum(At[j,]*log1p(X[j,]))+ sum(Bt[j,]*temp[(order2+j-1):j])+ sum(Ct[j,]*xreg[j,])
      temp   <- c(temp, vt[j])
      Int[j] <- mut[j,]
      AR[j] <- sum(At[j, ]*log1p(X[j, ]))
      MA[j] <- sum(Bt[j, ]*temp[(order2+j-1):j])
      for (k in 1:n){
        XREG[j,k] <- array(Ct[j,k] * xreg[j,k])
      }
    }
    
    comp <- data.frame(Int, AR, MA, XREG)
    
    data_long <- comp %>%
      dplyr::select(!c(Int))  %>%
      summarise_all(mean) %>%
      gather(variable, value) 
    
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
    
    xreg <- xreg[-(1:order),]
    n <- ncol(xreg)
    
    mut <- as.data.frame(coeff_list[["mut"]][["Mean"]])
    At <- as.data.frame(coeff_list[["At"]][["Mean"]])
    
    Ct_mat <- coeff_list[["Ct"]]
    Ct = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
    for(k in 1:n){
      Ct[,k] <- Ct_mat[[paste0("Ct", k)]][["Mean"]]
    }
    # Ct <- as.data.frame(Ct)
    
    AR = array(At * log1p(X))
    XREG = matrix(0, nrow = length(data)-order, ncol = ncol(xreg))
    for (k in 1:n){
      XREG[,k] <- array(Ct[,k] * xreg[,k])
    }
    colnames(XREG) <- colnames(xreg)
    
    comp2 = (array(At * log1p(X)))
    comp3 <- array(rowSums(Ct * xreg))
    vt_mat <- array(unlist(mut)) + array(unlist(comp2)) + array(comp3)
    
    comp <- data.frame(mut, AR, XREG)
    
    data_long <- comp %>%
      dplyr::select(!c("coeff_list...mut......Mean...")) %>%
      summarise_all(mean) %>%
      gather(variable, value) 
    
  }
  
  return(plot.contribution.summary(data_long))
}

plot.contribution.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long, aes(x = variable, y = value))+
    coord_flip() + 
    # sina plot: 
    geom_bar(stat="identity", aes(fill = value)) +
    # print the mean absolute value: 
    geom_text(data = unique(data_long[, c("variable", "value")]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) + 
    scale_fill_gradient(low="#FFCC33", high="#6600CC", 
                        breaks=c(0.8,3.5), labels = NULL) + #labels=c("Negative","Positive")
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), # remove axis line
          legend.position="bottom") + 
    # geom_hline(yintercept = 0) + # the vertical line
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "Posterior contribution (average impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}
