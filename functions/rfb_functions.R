############ Generic functions for neat MP functions
sample_index <- function(x, Data, reps, nyrs = 5) {
  Ind.ts <- Data@Ind[x, (length(Data@Ind[x, ])-nyrs+1):length(Data@Ind[x, ])]
  Ind5 <- trlnorm(reps * nyrs, Ind.ts, Data@CV_Ind[x])
  matrix(Ind5, ncol = reps)
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

# r = 1 for category 4 stocks
r_category4 <- function() return(1)





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
  Lmean <- numeric(length = dim(CAL)[2])
  Lc.ind <- which(LFS <= CAL_bins)[1]
  # Loop across years
  for(i in 1:length(Lmean)) {
    Lmean[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)] * CAL_bins[Lc.ind:length(CAL_bins)])/
      sum(CAL[i, Lc.ind:length(CAL_bins)])
  }
  # Stochastic reference points
  F01 <- numeric(length = length(Linf))
  Zcurr <- numeric(length = length(Linf))
  for(i in 1:length(Linf)) {
    F01[i] <- get_F01(Linf = Linf[i], K = K[i], Lc = LFS, 
                      M = M[i], wla = wla, wlb = wlb)
    Zcurr[i] <- model_select_GH(mlen = Lmean, K = K[i], Linf = Linf[i], Lc = LFS)
  }
    
  res <- F01/(Zcurr - M)
  res[is.infinite(res) | res < 0 | is.na(res)] <- 1
  return(res)
}

# f = F0.1/(Z-M) from GH (w/effort?)
f_GHeffort <- function(CAL, CAL_bins, effort, Linf, K, M, t0, wla, wlb, LFS) {
  Lmean <- numeric(length = dim(CAL)[2])
  Lc.ind <- which(LFS <= CAL_bins)[1]
  # Loop across years
  for(i in 1:length(Lmean)) {
    Lmean[i] <- sum(CAL[i, Lc.ind:length(CAL_bins)] * CAL_bins[Lc.ind:length(CAL_bins)])/
      sum(CAL[i, Lc.ind:length(CAL_bins)])
  }
  # Stochastic reference points
  F01 <- numeric(length = length(Linf))
  Zcurr <- numeric(length = length(Linf))
  for(i in 1:length(Linf)) {
    F01[i] <- get_F01(Linf = Linf[i], K = K[i], Lc = LFS, 
                      M = M[i], wla = wla, wlb = wlb)
    opt <- try(optim(c(1e-2, M[i]), GHeffort, method = "BFGS", 
                     Lbar = Lmean, ss = rep(1, length(Lmean)), 
                     eff = effort, Linf = Linf[i], K = K[i], a0 = t0[i], 
                     Lc = LFS, eff_init = effort[1], n_age = MaxAge,
                     control = list(maxit = 1e+06), hessian = FALSE), silent = TRUE)
    if(inherits(opt, "try-error")) Zcurr[i] <- NA
    else Zcurr[i] <- opt$par[1] * effort[length(effort)] + opt$par[2]
  }
  
  res <- F01/(Zcurr - M)
  res[is.infinite(res) | res < 0 | is.na(res)] <- 1
  return(res)
}


# f = spr/spr_target 
f_LBSPR <- function() {
}

# f = F0.1/(Z-M) from GH (w/effort?)



################ b-function 
# Iref (reference index) is the index at BMSY, calculated in the operating model.
# Thus, in ICES context, Ilim = 0.5 * Iref is the index at 0.5 BMSY.
# b = min(1, Icurrent/Itrigger) for category 3 stocks
# b = 1 for category 4, update once every 4 years
b_cat3 <- function(Icurr, Iref, w = 1.4) {
  Ilim <- 0.5 * Iref
  Itrigger <- w * Ilim
  Iratio <- Icurr/Itrigger
  b <- vapply(Iratio, function(x) min(x, 1), c(1))
  return(b)
}

b_cat4 <- function() return(0.8)