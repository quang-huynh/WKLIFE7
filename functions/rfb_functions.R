############ Generic functions for neat MP functions
sample_index <- function(x, Data, reps) {
  Ind.ts <- Data@Ind[x, ]
  Ind <- trlnorm(reps * length(Ind.ts), Ind.ts, Data@CV_Ind[x])
  matrix(Ind, ncol = reps)
}



############ r-functions

# r = exp(w * slope(log(I))), slope over most recent 5 years
calculate_slope <- function(log.ind) {
  dataframe <- data.frame(Ind = log.ind, Year = 1:length(log.ind))
  linreg <- lm(Ind ~ Year, data = dataframe)
  linreg$coefficients[2]
}

r_expIslope <- function(Ind5, w = 1) {
  slopes <- apply(log(Ind5), 2, calculate_slope)
  exp(w * slopes)
}

# r = 2 over 3 rule
r_2over3 <- function(Ind5, uncertainty_cap = FALSE, lower = 0.8, upper = 1.2) {
  I.num <- Ind5[4:5, ]
  I.den <- Ind5[1:3, ]
  if(ncol(Ind5) == 1) r <- mean(I.num, na.rm = TRUE)/mean(I.den, na.rm = TRUE)
  else r <- apply(I.num, 2, mean, na.rm = TRUE)/apply(I.den, 2, mean, na.rm = TRUE)
  if(uncertainty_cap) {
    r[r >= upper] <- upper
    r[r <= lower] <- lower
  }
  return(r)
}



############# f 
# f = Lmean/LF=M (Jardim et al. 2015) where Lc is the half modal length
# Lmean and Lc is deterministic, Linf/K/M/LF=M are stochastic
f_LBI <- function(CAL, CAL_bins, Linf, K, M) {
  modal.ind <- which.max(CAL)
  Lc.ind <- which(CAL[1:modal.ind] > 0.5 * CAL[modal.ind])[1]
  Lc <- CAL_bins[Lc.ind]
  Lmean <- sum(CAL[Lc.ind:length(CAL)] * CAL_bins[Lc.ind:length(CAL)])/sum(CAL[Lc.ind:length(CAL)])
  
  MK <- M/K
  LFM <- (Linf + 2 * MK * Lc)/(1 + 2 * MK)
  return(Lmean/LFM)
}

# f = M/(Z-M) from Beverton-Holt equation where Lc is the half modal length.
# Lmean and Lc is deterministic, Linf/K/M/LF=M are stochastic
f_BHE <- function(CAL, CAL_bins, Linf, K, M) {
  modal.ind <- which.max(CAL)
  Lc.ind <- which(CAL[1:modal.ind] > 0.5 * CAL[modal.ind])[1]
  Lc <- CAL_bins[Lc.ind]
  Lmean <- sum(CAL[Lc.ind:length(CAL)] * CAL_bins[Lc.ind:length(CAL)])/sum(CAL[Lc.ind:length(CAL)])
  
  Zcurr <- K * (Linf - Lmean)/(Lmean - Lc)
  res <- M/(Zcurr - M)
  res[is.infinite(res) | res < 0 | is.na(res)] <- 1
  return(res)
}

# f = F0.1/(Z-M) from GH (w/effort?)
f_GH <- function(CAL, CAL_bins, LFS, Linf, K, M, wla, wlb) {
  if(LFS > 0.95 * min(Linf)) LFS <- 0.95 * min(Linf)
  Lmean <- numeric(length = dim(CAL)[1])
  ss <- numeric(length = dim(CAL)[1])
  Lc.ind <- which(LFS <= CAL_bins)[1]
  # Loop across years
  for(i in 1:length(Lmean)) {
    Lmean[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)] * CAL_bins[Lc.ind:length(CAL_bins)])/
      sum(CAL[i, Lc.ind:length(CAL_bins)])
    ss[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)])
  }
  # Stochastic reference points
  F01 <- numeric(length = length(Linf))
  Zcurr <- numeric(length = length(Linf))
  for(i in 1:length(Linf)) {
    F01[i] <- get_F01(Linf = Linf[i], K = K[i], Lc = LFS, 
                      M = M[i], wla = wla, wlb = wlb)
    Zcurr[i] <- model_select_GH(mlen = Lmean, ss = ss, K = K[i], Linf = Linf[i], Lc = LFS)
  }
  res <- F01/(Zcurr - M)
  res[is.infinite(res) | res < 0 | is.na(res)] <- 1
  return(res)
}

# f = F0.1/F from GH
f_GHeffort <- function(CAL, CAL_bins, effort, Linf, K, M, t0, wla, wlb, LFS, MaxAge) {
  if(LFS > 0.95 * min(Linf)) LFS <- 0.95 * min(Linf)
  Lmean <- numeric(length = dim(CAL)[1])
  ss <- numeric(length = dim(CAL)[1])
  Lc.ind <- which(LFS <= CAL_bins)[1]
  # Loop across years
  for(i in 1:length(Lmean)) {
    Lmean[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)] * CAL_bins[Lc.ind:length(CAL_bins)])/
      sum(CAL[i, Lc.ind:length(CAL_bins)])
    ss[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)])
  }
  # Stochastic reference points
  F01 <- numeric(length = length(Linf))
  Fcurr <- numeric(length = length(Linf))
  effort <- effort/mean(effort, na.rm = TRUE)
  
  for(i in 1:length(Linf)) {
    opt <- try(optim(c(1e-2, M[i]), GHeffort, 
                     Lbar = Lmean, ss = rep(1, length(Lmean)), 
                     eff = effort, Linf = Linf[i], K = K[i], a0 = t0[i], 
                     Lc = LFS, eff_init = effort[1], n_age = MaxAge,
                     method = "BFGS", 
                     control = list(maxit = 1e+06), hessian = FALSE), silent = TRUE)
    #opt <- try(optim(1e-2, GHeffort_fixM, M = M[i], 
    #                 Lbar = Lmean, ss = ss, 
    #                 eff = effort, Linf = Linf[i], K = K[i], a0 = t0[i], 
    #                 Lc = LFS, eff_init = effort[1], n_age = MaxAge,
    #                 method = "BFGS",
    #                 control = list(maxit = 1e+06), hessian = FALSE), silent = TRUE)
    if(inherits(opt, "try-error") || opt$par[2] <= 0) {
      Fcurr[i] <- F01[i] <- NA
    } else {
      Fcurr[i] <- opt$par[1] * effort[length(effort)]
      F01[i] <- get_F01(Linf = Linf[i], K = K[i], Lc = LFS, M = opt$par[2], wla = wla, wlb = wlb)
    }
  }
  res <- F01/Fcurr
  res[is.infinite(res) | res < 0 | is.na(res)] <- 1
  return(res)
}


# f = spr/spr_target 
f_LBSPR <- function() {
  invisible()
}




################ b-function 
# Currently set Itrigger = minimum index in historical time series
b_cat3 <- function(Imatrix, w = 1.4) {
  Icurrent <- Imatrix[nrow(Imatrix), ]
  Ilim <- apply(Imatrix, 2, min, na.rm = TRUE)
  Itrigger <- w * Ilim
  Iratio <- Icurrent/Itrigger
  b <- vapply(Iratio, function(x) min(x, 1), c(1))
  return(b)
}


b_cat4 <- function(Data, buffer.interval = 4) {
  current.year <- max(Data@Year)
  LHYear <- Data@LHYear
  projection.year <- current.year - LHYear
  if(projection.year %% buffer.interval == 0) b <- 0.8 else b <- 1
  return(b)
}
