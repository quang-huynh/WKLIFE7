######### Category 4 HCR

######### LBI MPs
fLBI_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@CAL, Data@CAL_bins,
  Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort"
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_LBI(CAL = Data@CAL[x, dim(Data@CAL)[2], ], CAL_bins = CAL_bins, 
             Linf = Linfvec, K = Kvec, M = Mvec)
  
  b <- b_cat4(Data) 
  
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fLBI_cat4) <- "Output"
environment(fLBI_cat4) <- asNamespace("DLMtool")


############# BHE MPs

fBHE_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@CAL, Data@CAL_bins,
  Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort"
  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_BHE(CAL = Data@CAL[x, dim(Data@CAL)[2], ], CAL_bins = CAL_bins, 
            Linf = Linfvec, K = Kvec, M = Mvec)
  
  b <- b_cat4(Data) 
    
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fBHE_cat4) <- "Output"
environment(fBHE_cat4) <- asNamespace("DLMtool")


######### GH MPs
fGH_cat4 <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@CAL, Data@CAL_bins,
  Data@LFS, Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort"

  r <- 1
  
  Linfvec <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kvec <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  Mvec <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
  
  CAL_bins <- 0.5 * (Data@CAL_bins[1:(length(Data@CAL_bins) - 1)] + 
                       Data@CAL_bins[2:length(Data@CAL_bins)])
  f <- f_GH(CAL = Data@CAL[x, , ], CAL_bins = CAL_bins, LFS = Data@LFS[x], 
            Linf = Linfvec, K = Kvec, M = Mvec, wla = Data@wla[x], wlb = Data@wlb[x])
  
  b <- b_cat4(Data)
  
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(fGH_cat4) <- "Output"
environment(fGH_cat4) <- asNamespace("DLMtool")

