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
  
  if(!is.null(bsplines)){
    for (i in 1:ncol(summary)){
      mut[,i] <- bsplines %*% summary[,i]
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
        deltaA <- (summary_R[,i]+summary_K[,i])/2
        for(j in 1:order1){
          At[, j] <- bsplines_list[["bsplines"]] %*% deltaA[(j-1)*J+1:J]
        }
        At_mat[,i] <- At
        
        Bt <- matrix(0, length(data) - order, order2)
        deltaB <- (summary_R[,i]-summary_K[,i])/2
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
        deltaA <- (summary_R[,i]+summary_K[,i])/2
        for(j in 1:order1){
          At[, j] <- deltaA[j]
        }
        At_mat[,i] <- At
        
        Bt <- matrix(0, length(data) - order, order2)
        deltaB <- (summary_R[,i]-summary_K[,i])/2
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

predict <- function(data, xreg = NULL, coeff_list,
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
          vt[j,] <- mut_s[j] + sum(At_s[j]*log1p(X[j,]))+ sum(Bt_s[j]*temp_s[(order2+j-1):j])+ sum(Ct_s[j]*xreg[j,])
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

  return(data.frame(exp(vt_mat)))
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
          
          X <- NULL
          if(order > 0){
            for(j in order){
              ind <- (order-j+1):(length(data_pred) - j)
              X <- cbind(X, array(data_pred[ind]))
            } 
          }
          
          mut_k = mut_s %>%
            slice(rep(n(), each = k))
          mut_s = rbind(mut_s, mut_k)
          
          At_k = At_s %>% 
            slice(rep(n(), each = k))
          At_s = rbind(At_s,At_k)
          
          Bt_k = Bt_s %>% 
            slice(rep(n(), each = k))
          Bt_s = rbind(Bt_s,Bt_k)
          
          vt <- matrix(rep(0, (length(data_pred)-order)))
          for(j in 1:(length(data_pred)-order)){
            vt[j,] <- mut_s[j,] + sum(At_s[j,]*log1p(X[j, ]))+ sum(Bt_s[j,]*temp_s[(order2+j-1):j])
            temp_s    <- c(temp_s, vt[j,])
          } 
        
        t = length(temp_s)
        pred[k,i] <- as.numeric(exp(mut_s[t,] + sum(At_s[t,]*log1p(data[t]))+ sum(Bt_s[t,]*temp_s[(t)])))
        
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
            slice(rep(n(), each = k))
          mut_s = rbind(mut_s, mut_k)
          
          At_k = At_s %>% 
            slice(rep(n(), each = k))
          At_s = rbind(At_s,At_k)
          
          Bt_k = Bt_s %>% 
            slice(rep(n(), each = k))
          Bt_s = rbind(Bt_s,Bt_k)
          
          Ct_k = Ct_s %>% 
            slice(rep(n(), each = k))
          Ct_s = rbind(Ct_s,Ct_k)
          
          vt <- matrix(rep(0, (length(data_pred)-order)))
          for(j in 1:(length(data_pred)-order)){
            vt[j,] <- mut_s[j,] + sum(At_s[j, ]*log1p(X[j, ]))+ sum(Bt_s[j, ]*temp_s[(order2+j-1):j]) + sum(Ct_s[j,]*xreg_pred[j,])
            temp_s    <- c(temp_s, vt[j,])
          } 
        
        t = length(temp_s)
        pred[k,i] <- as.numeric(exp(mut_s[t,] + sum(At_s[t, ]*log1p(data[t]))+ sum(Bt_s[t, ]*temp_s[(t)]) + sum(Ct_s[t,]*xreg_pred[t,])))
        
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
        Ct_k = data.frame(Ct) %>% filter(row_number()==n()) %>% pull()
        data_pred <- data.frame(data) %>% filter(row_number()==n()) %>% pull()
        xreg_pred <- data.frame(xreg) %>% filter(row_number()==n()) %>% pull()
        for (k in 1:nrow(pred)){
        
          pred[k,i] <- as.numeric(exp(array(mut_k) + array(At_k* log1p(data_pred[k])) + array(Ct_k * xreg_pred)))
          data_pred <- c(data_pred, pred[k,i])
          
        }
      }
      p <- data.frame(forecast_date = seq(as.Date(lastDate)+1, by = "day", length.out = n.ahead))
      pred <- cbind(p, pred)
    }
  }
  return(pred)
}


  