######### Yield per recruit functions
# Get the yield per recruit as a function of F
get_YPR <- function(F_mort, Linf, K, Lc, M, wla, wlb) {
  Z <- F_mort + M
  Winf <- wla * Linf ^ wlb
  biomass <- (Winf/K) * beta(wlb+1, Z/K) * (1 - pbeta(Lc/Linf, wlb+1, Z/K)) * 
    (1 - Lc/Linf)^(Z/K)
  return(F_mort * biomass)
}

# Calculate the derivative of the YPR with respect to F, i.e. dY/dF
get_deriv <- function(F_mort, Linf, K, Lc, M, wla, wlb) {
  numDeriv::grad(get_YPR, x = F_mort, Linf = Linf, K = K, Lc = Lc, M = M, 
                 wla = wla, wlb = wlb)
}

F01_rootfinder <- function(F_mort, Linf, K, Lc, M, wla, wlb) {
  slope.origin <- get_deriv(F_mort = 0, Linf = Linf, K = K, Lc = Lc, M = M, 
                            wla = wla, wlb = wlb)
  slope.F <- get_deriv(F_mort = F_mort, Linf = Linf, K = K, Lc = Lc, M = M, 
                       wla = wla, wlb = wlb)
  slope.F - 0.1 * slope.origin
}

get_F01 <- function(Linf, K, Lc, M, wla, wlb, maxF = 3) {
  res <- try(uniroot(F01_rootfinder, interval = c(1e-3, maxF), Linf = Linf, K = K, Lc = Lc, M = M, 
          wla = wla, wlb = wlb)$root)
  if(inherits(res, 'try-error')) res <- NA
  return(res)
}

# Loops through 0 - 3 changes in mortality
model_select_GH <- function(mlen, ss = rep(1, length(mlen)), K, Linf, Lc) {
  mlen[mlen <= 0 | is.na(mlen)] <- -99
  ss[mlen == -99] <- 0
  
  #if(length(mlen) > 40) { # Use the most recent 40 years of data
  #  yrind <- (length(mlen)-39):length(mlen)
  #  mlen <- mlen[yrind]
  #  ss <- ss[yrind]
  #}
  for(i in 0:5) {
    if(i == 0) {
      stpar <- 0.5
      res0 <- try(optim(stpar, GHeq_LL, method = "BFGS", 
                        Lbar = mlen, ss = ss,
                        K = K, Linf = Linf, Lc = Lc, control = list(maxit = 1e+06), 
                        hessian = FALSE), silent = TRUE)
      if(inherits(res0, 'try-error'))
        aic0 <- NA
      else aic0 <- 2 * (res0$value + 2)
    }
    if(i > 0) {
      stpar <- c(rep(0.5, i+1), c(1:i) * length(mlen)/(i+1))
      res <- try(optim(stpar, GH_LL, method = "BFGS", 
                       Lbar = mlen, ss = ss, nbreaks = i, 
                       K = K, Linf = Linf, Lc = Lc, control = list(maxit = 1e+06), 
                       hessian = FALSE), silent = TRUE)
      if(inherits(res, 'try-error')) {
        assign(paste0('res', i), NA)
        assign(paste0('aic', i), NA)
      } else {
        npar <- 2 * i + 2
        assign(paste0('res', i), res)
        assign(paste0('aic', i), 2 * getElement(get(paste0('res', i)), 'value') + npar)
      }
    }
  }
  
  testvec <- c(aic0, aic1, aic2, aic3, aic4, aic5)
  if(all(is.na(testvec))) 
    Z.ans <- NA
  else {
    index.best.nbreak <- which.min(testvec) # Find the model with the minimum AIC
    if(index.best.nbreak > 1) {
      # Choose model with n-1 change points instead of n change points if the difference in AIC is less than 2.
      choice1 <- testvec[index.best.nbreak] + 2
      choice2 <- testvec[index.best.nbreak - 1]
      if(!is.na(choice2)) {
        if(choice1 > choice2) index.best.nbreak <- index.best.nbreak-1 
      }
    }
    best.model <- get(paste0('res', index.best.nbreak - 1))
    Z.ans <- best.model$par[index.best.nbreak]
    if(Z.ans < 0) Z.ans <- NA
  }
  
  return(Z.ans)
}

