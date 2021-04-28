#' @title The functions to fit time-constant AR(p), INGARCH(p,q), 
#'      ARX(p), INGARCHX(p,q) model for count data
#' @description Takes the count-valued time series as input and returns MCMC samples
#' @references 
#'
#' @param data is the time-series of count-valued data
#' @param order is the order of the AR component for AR(p) model
#' @param order1 is the order of the AR component for INGARCH(p,q) model
#' @param order2 is the order of the MA component for INGARCH(p,q) model
#' @param Total_itr is the total number of iterations of MCMC
#' @param burn is the number of burn-in MCMC samples

#' @return fit._.log returns a list of the posterior samples of parameters
#' (deltamufn, sigma2latfn, deltaRfn, deltaKfn, deltaCfn), predictions (vtp) 
#' and MAPE (pred) \cr

fit.AR.log <- function(data, order = 1, Total_itr = 5000, burn = 2500){
  set.seed(1)
  
  Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    vt <- array(mut) + array(comp2)
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    vt <- array(mut) + array(comp2) 
    
    compo   <- exp(vt) - Y 
    
    return(sum(array(compo))+ x/100)
  }
  
  UR <- function(x){
    xc <- x
    deltaA <- (1-exp(xc))/(1+exp(xc))
    
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    vt <- array(mut) + array(comp2) 
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(xc^2) / (2*100))
  }
  
  grad_UR <- function(x, tempsig=sigma2lat){
    xc <- x
    deltaA <- (1-exp(xc))/(1+exp(xc))
    
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At *log1p(X)))}
    
    vt <- array(mut) + array(comp2) 
    
    if(order==1){compo   <- array(array((exp(vt) - Y)*(array(log1p(X)))))}
    if(order > 1){
      compo <- NULL
      for(j in 1:order){
        compo   <- c(compo, array(array((exp(vt) - Y)*array(log1p(X[, j])))))
      }
    }
    
    xcder <- - exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc)) 
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  # Data and priors 
  
  Y <- data
  if(order>0){Y <- data[-(1:order)]}
  X <- NULL
  if(order > 0){
    for(j in 1:order){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, array(data[ind]))
    } 
  }
  
  At <- 0
  fit.G <- tsglm(Y, model=list(past_obs=1), distr = "poisson", link = "log")
  deltamu <- fit.G$coefficients[1]
  deltaA <- fit.G$coefficients[2]

  deltaR <- log((1-deltaA)/(1+deltaA))
  
  # At <- 0
  # deltamu <- rnorm(1)
  # deltaA <- runif(1,-1,1)
  # 
  # deltaR <- log((1-deltaA)/(1+deltaA))
  
  mut     <- rep(deltamu, (length(data)-order))
  
  if(order>0){
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
  }
  
  comp2 = array(At * log1p(X))
  if(order==0){comp2 = rep(0, length(mut))}
  if(order>1){comp2 = array(rowSums(At * log1p(X)))}
  
  vt <- array(mut) + array(comp2) 
  
  # Set initial values and run HMC sampling
  
  itr <- 0
  armu <- 0
  sdmu <-  1e-4
  Lmu <- 1
  arR <- 0
  sdR <- 1e-4
  LR <-  1
  vts <- list()
  deltamuls <- list()
  deltaRls <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  
  while(itr < Total_itr){
    
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, Lmu, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut     <- rep(deltamu, (length(data)-order))
    
    temp   <- HMC(UR, grad_UR, sdR, LR, deltaR, arR)
    deltaR <- temp$up
    arR <- temp$arc
    
    deltaA <- (1-exp(deltaR))/(1+exp(deltaR))
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    
    vt <- array(mut) + array(comp2) 
    
    Pred[itr] <- MLmetrics::MAPE(exp(vt), Y)
    vts[[itr]] <- vt
    
    if(itr %% 100 == 0){
      if(order > 0){
        
        ar <- arR/ itr
        cat(ar, "acceptance rate for R")
        if(ar<.45){sdR <- sdR * (.1)}
        if(ar>.70){sdR <- sdR * (10)}
        
      }
      
      
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.35){sdmu <- sdmu * (.1)}
      if(ar>.70){sdmu <- sdmu * (10)}
    }
    
    deltamuls[[itr]] <- deltamu
    deltaRls[[itr]] <- deltaR
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    plot(exp(vt))
  }
  close(pb)
  out <- list(vtp = vts[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr],
              deltamufn = deltamuls[(burn+1):Total_itr], 
              deltaRfn = deltaRls[(burn+1):Total_itr])
  return(out)
}

fit.INGARCH.log <- function(data, order1 = 1, order2 = 1, Total_itr = 5000, burn = 2500){
  set.seed(1)
  order = max(order1,order2)
  
  Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    sigma2m <- Bt
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y 
    
    return(sum(array(compo))+ x/100)
  }
  
  UR <- function(x){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (x+deltaK)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (x-deltaK)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(xc^2) / (2*100))
  }
  
  grad_UR <- function(x, deltaKk=deltaK, tempsig=sigma2lat){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (x+deltaKk)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (x-deltaKk)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At *log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- tempsig #sigma2lat
      sigma2m <- Bt
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vt[i])
      } 
    }
    
    if(order1==1){compo   <- array(1/2 * array((exp(vt) - Y)*(array(log1p(X)+sigma2m))))}
    if(order1 > 1){
      compo <- NULL
      for(j in 1:order1){
        compo   <- c(compo, array(1/2 * array((exp(vt) - Y)*array(log1p(X[, j])+sigma2m[, j]))))
      }
    }
    xcder <- -exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc)) 
    
    #expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    #return(array(expxdev %*% compo) + (x) / (100))
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  UK <- function(x){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (deltaR+x)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (deltaR-x)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo) + sum(xc^2) / (2*100))
  }
  
  grad_UK <- function(x, deltaRr = deltaR, tempsig=sigma2lat){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (deltaRr+x)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (deltaRr-x)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    sigma2m <- Bt
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vt[i])
      } 
    }
    
    if(order2==1){compo   <- array(1/2 * array((exp(vt)-Y)*array(log1p(X)-sigma2m)))}
    if(order2 > 1){
      compo <- NULL
      for(j in 1:order2){
        compo   <- c(compo, array(1/2 * array((exp(vt)-Y)*array(log1p(X[,j])-sigma2m[, j]))))
      }
    }
    xcder <- -exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc))  
    #expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    #return(array(expxdev %*% compo) + (x) / (100))
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  
  US <- function(x, deltaRr = deltaR, deltaKk = deltaK){
    deltaA <- (deltaRr+deltaKk)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j] 
    }
    deltaB <- (deltaRr-deltaKk)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    
    vt <- length(Y)
    temp <- x
    if(order2>0){
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
        #sigma2m[i, ] <- temp[(order2+i-1):i]
      } 
    }
    compo   <- sum(exp(vt) - Y * vt) + sum(exp(x) - (data[1:order])*x) #+ sum(x - (data[1:order])*log(x))
    return(sum(compo) + sum(x^2) / (2*100))
  }
  
  grad_US <- function(x, deltaRr, deltaKk){
    ret <- jacobian(US, x, deltaRr = deltaR, deltaKk = deltaK)
    return(array(ret))
  }
  
  # Data and priors
  
  Y <- data
  if(order>0){Y <- data[-(1:order)]}
  X <- NULL
  if(order1 > 0){
    for(j in 1:order1){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, array(data[ind]))
    } 
  }
  
  At <- 0
  
  fit.G <- tsglm(data, model=list(past_obs=1, past_mean=1), distr = "poisson", link = "log")
  deltamu <- fit.G$coefficients[1]
  deltaA <- fit.G$coefficients[2]
  deltaB <- fit.G$coefficients[3]
  
  deltaR <- (deltaA+deltaB)/2
  deltaK <- (deltaA-deltaB)/2
  
  deltaRu <- log((1-deltaR)/(1+deltaR))
  deltaKu <- log((1-deltaK)/(1+deltaK))
  
  mut     <- rep(deltamu, (length(data)-order))
  
  if(order1>0){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
  }
  
  if(order2>0){
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
  }
  
  sigma2lat <- rnorm(1)#
  
  delta <- c(sigma2lat, deltaRu, deltaKu)
  temp <- sigma2lat
  if(order2>0){
    vt <- rep(0, length(data)-order) 
    for(i in 1:(length(data)-order)){
      vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
      temp    <- c(temp, vt[i])
    } 
  }
  
  # Set initial values and run HMC sampling
  
  itr <- 0
  armu <- 0
  sdmu <-  1e-4 
  Lmu <- 1
  arS <- 0
  sdS <-  1e-4 
  LS <- 1
  arRu <- 0
  sdRu <-  1e-4 
  LRu <- 1
  arKu <- 0
  sdKu <- 1e-4
  LKu <-  1
  vts <- list()
  deltamuls <- list()
  sigma2latls <- list()
  deltaRls <- list()
  deltaKls <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  
  while(itr < Total_itr){
    
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, Lmu, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut     <- rep(deltamu, (length(data)-order))
    
    temp <- HMC(US, grad_US, sdS, LS, sigma2lat, arS)
    sigma2lat <- temp$up
    arS <- temp$arc
    
    temp <- HMC(UR, grad_UR, sdRu, LRu, deltaRu, arRu)
    deltaRu <- temp$up
    arRu <- temp$arc
    deltaR <- (1-exp(deltaRu))/(1+exp(deltaRu))
    
    temp <- HMC(UK, grad_UK, sdKu, LKu, deltaKu, arKu)
    deltaKu <- temp$up
    arKu <- temp$arc
    deltaK <- (1-exp(deltaKu))/(1+exp(deltaKu))
    
    deltaA <- (deltaR+deltaK)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    deltaB <- (deltaR-deltaK)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i]) 
        temp    <- c(temp, vt[i])
      } 
    }
   
    Pred[itr] <- MLmetrics::MAPE(exp(vt), Y)
    vts[[itr]] <- vt
    
    if(itr %% 100 == 0){
      if(order > 0){
        
        ar <- arS/ itr
        cat(ar, "acceptance rate for S")
        if(ar<.45){sdS <- sdS * (.5)}
        if(ar>.70){sdS <- sdS * (5)}
        
        ar <- arRu/ itr
        cat(ar, "acceptance rate for R")
        if(ar<.45){sdRu <- sdRu * (.5)}
        if(ar>.70){sdRu <- sdRu * (5)}
        
        ar <- arKu/ itr
        cat(ar, "acceptance rate for K")
        if(ar<.45){sdKu <- sdKu * (.5)}
        if(ar>.70){sdKu <- sdKu * (5)}
        
      }
      
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.45){sdmu <- sdmu * (.5)}
      if(ar>.70){sdmu <- sdmu * (5)}
    }
    
    deltamuls[[itr]] <- deltamu
    sigma2latls[[itr]] <- sigma2lat
    deltaRls[[itr]] <- deltaR
    deltaKls[[itr]] <- deltaK
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    plot(exp(vt))
  }
  close(pb)
  out <- list(vtp = vts[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr],
              deltamufn = deltamuls[(burn+1):Total_itr], sigma2latfn = sigma2latls[(burn+1):Total_itr], 
              deltaRfn = deltaRls[(burn+1):Total_itr], deltaKfn = deltaKls[(burn+1):Total_itr])
  return(out)
}

fit.ARX.log <- function(data, Xreg, order = 1, Total_itr = 5000, burn = 2500){
  set.seed(1)
  
  Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    compo   <- exp(vt) - Y 
    
    return(sum(array(compo))+ x/100)
  }
  
  UR <- function(x){
    xc <- x
    deltaA <- (1-exp(xc))/(1+exp(xc))
    
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(xc^2) / (2*100))
  }
  
  grad_UR <- function(x){
    xc <- x
    deltaA <- (1-exp(xc))/(1+exp(xc))
    
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At *log1p(X)))}
    
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    if(order==1){compo   <- array(array((exp(vt) - Y)*(array(log1p(X)))))}
    if(order > 1){
      compo <- NULL
      for(j in 1:order){
        compo   <- c(compo, array(array((exp(vt) - Y)*array(log1p(X[, j])))))
      }
    }
    
    xcder <- - exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc)) 
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  UC <- function(x){
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <- x[j]
    }
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_UC <- function(x, deltaAa=deltaA){
    
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaAa[j]
    }
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <-x[j]
    }
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    compo <- NULL
    for(j in 1:n){
      compo   <- c(compo, array(array((Xreg[, j]*exp(vt))-(Y*Xreg[, j]))))
    }
    return(array(compo)+ x / 100)
  }
  
  # Data and priors 
  
  Y <- data
  n <- ncol(Xreg)
  if(order>0){Y <- data[-(1:order)]
  Xreg <- Xreg[-(1:order),]
  }
  X <- NULL
  if(order > 0){
    for(j in 1:order){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, array(data[ind]))
    } 
  }
  
  At <- 0
  fit.G <- tsglm(Y, model=list(past_obs=1), distr = "poisson", link = "log", xreg = Xreg)
  deltamu <- fit.G$coefficients[1]
  deltaA <- fit.G$coefficients[2]
  C = c(); deltaC = c()
  for (i in 1:n){
    C[i] <- fit.G$coefficients[i+2]
    deltaC <- c(deltaC, C[i])
  }
  
  deltaR <- log((1-deltaA)/(1+deltaA))
  mut     <- rep(deltamu, (length(data)-order))
  
  if(order>0){
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
  }
  
  Ct <- matrix(0, length(data) - order, n)
  for(j in 1:n){
    Ct[, j] <- deltaC[j]
  }
  
  comp2 = array(At * log1p(X))
  if(order==0){comp2 = rep(0, length(mut))}
  if(order>1){comp2 = array(rowSums(At * log1p(X)))}
  
  comp3 <- array(rowSums(Ct * Xreg))
  vt <- array(mut) + array(comp2) + array(comp3)
  
  # Set initial values and run HMC sampling
  
  itr <- 0
  armu <- 0
  sdmu <-  1e-4
  Lmu <- 1
  arR <- 0
  sdR <- 1e-4
  LR <-  1
  arC <- 0
  sdC <- 1e-4 
  LC <- 1
  vts <- list()
  deltamuls <- list()
  deltaRls <- list()
  deltaCls <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  
  while(itr < Total_itr){
    
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, Lmu, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut     <- rep(deltamu, (length(data)-order))
    
    temp   <- HMC(UR, grad_UR, sdR, LR, deltaR, arR)
    deltaR <- temp$up
    arR <- temp$arc
    
    deltaA <- (1-exp(deltaR))/(1+exp(deltaR))
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- deltaA[j]
    }
    comp2 = array(At * log1p(X))
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * log1p(X)))}
    
    temp    <- HMC(UC, grad_UC, sdC, LC, deltaC, arC)
    deltaC <- temp$up
    arC     <- temp$arc
    
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <- deltaC[j]
    }
    
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    Pred[itr] <- MLmetrics::MAPE(exp(vt), Y)
    vts[[itr]] <- vt
    
    if(itr %% 100 == 0){
      if(order > 0){
        
        ar <- arR/ itr
        cat(ar, "acceptance rate for combo")
        if(ar<.45){sdR <- sdR * (.1)}
        if(ar>.70){sdR <- sdR * (10)}
        
      }
      
      ar <- arC/ itr
      cat(ar, "acceptance rate for C")
      if(ar<.45){sdC <- sdC * (.1)}
      if(ar>.70){sdC <- sdC * (10)}
      
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.35){sdmu <- sdmu * (.1)}
      if(ar>.70){sdmu <- sdmu * (10)}
    }
    
    deltamuls[[itr]] <- deltamu
    deltaRls[[itr]] <- deltaR
    deltaCls[[itr]] <- deltaC
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    plot(exp(vt))
  }
  close(pb)
  out <- list(vtp = vts[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr],
              deltamufn = deltamuls[(burn+1):Total_itr], 
              deltaRfn = deltaRls[(burn+1):Total_itr], 
              deltaCfn = deltaCls[(burn+1):Total_itr])
  return(out)
}

fit.INGARCHX.log <- function(data, Xreg, order1 = 1, order2 = 1, Total_itr = 5000, burn = 2500){
  set.seed(1)
  order = max(order1,order2)
  
  Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    sigma2m <- Bt
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut <- array(rep(x, (length(data)-order)))
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y 
    
    return(sum(array(compo))+ x/100)
  }
  
  UR <- function(x){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (x+deltaK)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (x-deltaK)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(xc^2) / (2*100))
  }
  
  grad_UR <- function(x, deltaKk=deltaK, tempsig=sigma2lat){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (x+deltaKk)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (x-deltaKk)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At *log1p(X)))}
    
    # vt <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- tempsig #sigma2lat
      sigma2m <- Bt
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vt[i])
      } 
    }
    
    if(order1==1){compo   <- array(1/2 * array((exp(vt) - Y)*(array(log1p(X)+sigma2m))))}
    if(order1 > 1){
      compo <- NULL
      for(j in 1:order1){
        compo   <- c(compo, array(1/2 * array((exp(vt) - Y)*array(log1p(X[, j])+sigma2m[, j]))))
      }
    }
    xcder <- - exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc)) 
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  UK <- function(x){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (deltaR+x)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (deltaR-x)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo) + sum(xc^2) / (2*100))
  }
  
  grad_UK <- function(x, deltaRr = deltaR, tempsig=sigma2lat){
    xc <- x
    x <- (1-exp(xc))/(1+exp(xc))
    deltaA <- (deltaRr+x)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    
    deltaB <- (deltaRr-x)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    # vt <- array(mut) + array(comp2) 
    
    sigma2m <- Bt
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vt[i])
      } 
    }
    
    if(order2==1){compo   <- array(1/2 * array((exp(vt)-Y)*array(log1p(X)-sigma2m)))}
    if(order2 > 1){
      compo <- NULL
      for(j in 1:order2){
        compo   <- c(compo, array(1/2 * array((exp(vt)-Y)*array(log1p(X[,j])-sigma2m[, j]))))
      }
    }
    xcder <- - exp(xc)*(1-exp(xc))/(1+exp(xc))^2 - exp(xc)/(1+exp(xc)) 
    return(sum(array(compo)*xcder)+ xc/100)
  }
  
  UC <- function(x){
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <- x[j]
    }
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])+ sum(Ct[i,]*Xreg[i,])
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo   <- exp(vt) - Y * vt
    return(sum(compo)+ sum(x^2) / (2*100))
  }
  
  grad_UC <- function(x, deltaRr=deltaR, deltaKk=deltaK, tempsig=sigma2lat){
    
    deltaA <- (deltaRr+deltaKk)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    deltaB <- (deltaRr-deltaKk)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <-x[j]
    }
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    sigma2m <- Bt
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])+ sum(Ct[i,]*Xreg[i,])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vt[i])
      } 
    }
    
    compo <- NULL
    for(j in 1:n){
      compo   <- c(compo, array(array((Xreg[, j]*exp(vt))-(Y*Xreg[, j]))))
    }
    return(array(compo)+ x / 100)
  }
  
  
  US <- function(x, deltaRr = deltaR, deltaKk = deltaK){
    deltaA <- (deltaRr+deltaKk)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j] 
    }
    deltaB <- (deltaRr-deltaKk)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    
    vt <- length(Y)
    temp <- x
    if(order2>0){
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vt[i])
      } 
    }
    compo   <- sum(exp(vt) - Y * vt) + sum(exp(x) - (data[1:order])*x) #+ sum(x - (data[1:order])*log(x))
    return(sum(compo) + sum(x^2) / (2*100))
  }
  
  grad_US <- function(x, deltaRr, deltaKk){
    ret <- jacobian(US, x, deltaRr = deltaR, deltaKk = deltaK)
    return(array(ret))
  }
  
 # Data and priors
  
  Y <- data
  n <- ncol(Xreg)
  if(order>0){Y <- data[-(1:order)]
  Xreg <- Xreg[-(1:order),]
  }
  X <- NULL
  if(order > 0){
    for(j in 1:order1){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, array(data[ind]))
    } 
  }
  
  At <- 0
  fit.G <- tsglm(Y, model=list(past_obs=1, past_mean=1), distr = "poisson", link = "log", xreg = Xreg)
  deltamu <- fit.G$coefficients[1]
  deltaA <- fit.G$coefficients[2]
  deltaB <- fit.G$coefficients[3]
  C = c(); deltaC = c()
  for (i in 1:n){
    C[i] <- fit.G$coefficients[i+3]
    deltaC <- c(deltaC, C[i])
  }
  
  deltaR <- (deltaA+deltaB)/2
  deltaK <- (deltaA-deltaB)/2
  
  deltaRu <- log((1-deltaR)/(1+deltaR))
  deltaKu <- log((1-deltaK)/(1+deltaK))
  
  mut     <- rep(deltamu, (length(data)-order))
  
  if(order1>0){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
  }
  
  if(order2>0){
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
  }
  
  sigma2lat <- rnorm(1)#
  
  delta <- c(sigma2lat, deltaRu, deltaKu)
  temp <- sigma2lat
  if(order2>0){
    vt <- rep(0, length(data)-order) 
    for(i in 1:(length(data)-order)){
      vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i])
      temp    <- c(temp, vt[i])
    } 
  }
  
  # Set initial values and run HMC sampling
  
  itr <- 0
  armu <- 0
  sdmu <-  1e-4 
  Lmu <- 1
  arS <- 0
  sdS <-  1e-4 
  LS <- 1
  arRu <- 0
  sdRu <-  1e-4 
  LRu <- 1
  arKu <- 0
  sdKu <- 1e-4
  LKu <-  1
  arC <- 0
  sdC <- 1e-4
  LC <- 1
  vts <- list()
  deltamuls <- list()
  sigma2latls <- list()
  deltaRls <- list()
  deltaKls <- list()
  deltaCls <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  
  while(itr < Total_itr){
    
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, Lmu, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut     <- rep(deltamu, (length(data)-order))
    
    temp <- HMC(US, grad_US, sdS, LS, sigma2lat, arS)
    sigma2lat <- temp$up
    arS <- temp$arc
    
    temp <- HMC(UR, grad_UR, sdRu, LRu, deltaRu, arRu)
    deltaRu <- temp$up
    arRu <- temp$arc
    deltaR <- (1-exp(deltaRu))/(1+exp(deltaRu))
    
    temp <- HMC(UK, grad_UK, sdKu, LKu, deltaKu, arKu)
    deltaKu <- temp$up
    arKu <- temp$arc
    deltaK <- (1-exp(deltaKu))/(1+exp(deltaKu))
    
    deltaA <- (deltaR+deltaK)/2
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- deltaA[j]
    }
    deltaB <- (deltaR-deltaK)/2
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- deltaB[j]
    }
    
    comp2 = array(At * log1p(X))
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * log1p(X)))}
    
    temp    <- HMC(UC, grad_UC, sdC, LC, deltaC, arC)
    deltaC <- temp$up
    arC     <- temp$arc
    
    Ct <- matrix(0, length(data) - order, n)
    for(j in 1:n){
      Ct[, j] <- deltaC[j]
    }
    
    comp3 <- array(rowSums(Ct * Xreg))
    vt <- array(mut) + array(comp2) + array(comp3)
    
    # vt <- array(mut) + array(comp2) 
    
    if(order2>0){
      temp <- sigma2lat
      vt <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vt[i] <- mut[i] + sum(At[i, ]*log1p(X[i, ]))+ sum(Bt[i, ]*temp[(order2+i-1):i]) + sum(Ct[i,]*Xreg[i,])
        temp    <- c(temp, vt[i])
      } 
    }
    
    Pred[itr] <- MLmetrics::MAPE(exp(vt), Y)
    vts[[itr]] <- vt
    
    if(itr %% 100 == 0){
      if(order > 0){
        
        ar <- arS/ itr
        cat(ar, "acceptance rate for combo")
        if(ar<.45){sdS <- sdS * (.5)}
        if(ar>.70){sdS <- sdS * (5)}
        
        ar <- arRu/ itr
        cat(ar, "acceptance rate for combo")
        if(ar<.45){sdRu <- sdRu * (.5)}
        if(ar>.70){sdRu <- sdRu * (5)}
        
        ar <- arKu/ itr
        cat(ar, "acceptance rate for combo")
        if(ar<.45){sdKu <- sdKu * (.5)}
        if(ar>.70){sdKu <- sdKu * (5)}
        
      }
      
      ar <- arC/ itr
      cat(ar, "acceptance rate for C")
      if(ar<.45){sdC <- sdC * (.1)}
      if(ar>.70){sdC <- sdC * (10)}
      
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.35){sdmu <- sdmu * (.1)}
      if(ar>.70){sdmu <- sdmu * (10)}
    }
    
    deltamuls[[itr]] <- deltamu
    sigma2latls[[itr]] <- sigma2lat
    deltaRls[[itr]] <- deltaR
    deltaKls[[itr]] <- deltaK
    deltaCls[[itr]] <- deltaC
    
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    plot(exp(vt))
  
  }
  close(pb)
  out <- list(vtp = vts[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr],
              deltamufn = deltamuls[(burn+1):Total_itr], sigma2latfn = sigma2latls[(burn+1):Total_itr], 
              deltaRfn = deltaRls[(burn+1):Total_itr], deltaKfn = deltaKls[(burn+1):Total_itr],
              deltaCfn = deltaCls[(burn+1):Total_itr])
  return(out)
}
