######### Category 4 HCR

######### LBI MPs
fLBI_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind, Data@CAL, Data@CAL_bins,
  Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort,
  Data@Iref, Data@CV_Iref"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_LBI(CAL = Data@CAL[x, dim(Data@CAL)[2], ], CAL_bins = CAL_bins, 
             Linf = Linfvec, K = Kvec, M = Mvec)
  
  Icurr.vec <- Ind5[5, ]
  Iref.vec <- trlnorm(reps, Data@Iref[x], Data@CV_Iref[x])
  b <- b_cat4()
  
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fLBI_cat4) <- "Output"
environment(fLBI_cat4) <- asNamespace("DLMtool")


############# BHE MPs

fBHE_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind, Data@CAL, Data@CAL_bins,
  Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort,
  Data@Iref, Data@CV_Iref"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_BHE(CAL = Data@CAL[x, dim(Data@CAL)[2], ], CAL_bins = CAL_bins, 
            Linf = Linfvec, K = Kvec, M = Mvec)
  
  Icurr.vec <- Ind5[5, ]
  Iref.vec <- trlnorm(reps, Data@Iref[x], Data@CV_Iref[x])
  b <- b_cat4()
    
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fBHE_cat4) <- "Output"
environment(fBHE_cat4) <- asNamespace("DLMtool")


######### GH MPs
fGH_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort,
  Data@Iref, Data@CV_Iref"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_GH(CAL = Data@CAL[x, , ], CAL_bins = CAL_bins, LFS = Data@LFS[x], 
            Linf = Linfvec, K = Kvec, M = Mvec, wla = Data@wla[x], wlb = Data@wlb[x])
  
  Icurr.vec <- Ind5[5, ]
  Iref.vec <- trlnorm(reps, Data@Iref[x], Data@CV_Iref[x])
  b <- b_cat4()
    
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fGH_cat4) <- "Output"
environment(fGH_cat4) <- asNamespace("DLMtool")


# GH effort

fGHe_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort,
  Data@Iref, Data@CV_Iref"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  t0vec <- trlnorm(reps, abs(Data@vbt0[x]), Data@CV_vbt0[x])
  if(Data@vbt0[x] < 0) t0vec <- -1 * t0vec
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_GHeffort(CAL = Data@CAL[x, , ], CAL_bins = CAL_bins, LFS = Data@LFS[x], 
                  effort = Data@Cat[x, ]/Data@Ind[x, ],
                  Linf = Linfvec, K = Kvec, M = Mvec, t0 = t0vec, 
                  wla = Data@wla[x], wlb = Data@wlb[x])
  
  Icurr.vec <- Ind5[5, ]
  Iref.vec <- trlnorm(reps, Data@Iref[x], Data@CV_Iref[x])
  b <- b_cat4()
  
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fGHe_cat4) <- "Output"
environment(fGHe_cat4) <- asNamespace("DLMtool")


r5sl_fGHe_b3 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort,
  Data@Iref, Data@CV_Iref"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- r_expIslope(Ind5)
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  t0vec <- trlnorm(reps, abs(Data@vbt0[x]), Data@CV_vbt0[x])
  if(Data@vbt0[x] < 0) t0vec <- -1 * t0vec
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_GHeffort(CAL = Data@CAL[x, , ], CAL_bins = CAL_bins, LFS = Data@LFS[x], 
                  effort = Data@Cat[x, ]/Data@Ind[x, ],
                  Linf = Linfvec, K = Kvec, M = Mvec, t0 = t0vec, 
                  wla = Data@wla[x], wlb = Data@wlb[x])
  
  Icurr.vec <- Ind5[5, ]
  Iref.vec <- trlnorm(reps, Data@Iref[x], Data@CV_Iref[x])
  b <- b_cat3(Icurr.vec, Iref.vec)
  
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(r5sl_fGHe_b3) <- "Output"
environment(r5sl_fGHe_b3) <- asNamespace("DLMtool")
