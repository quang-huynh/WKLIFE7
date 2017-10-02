# Use the 'source()' function in order to load all functions into the workspace.

##### A proposal for nomenclature of control rules.
#
# For the control rule of the format: TAC = Catch * r * f * b, the function name of
# the control rule is: x_y_z where
#
# Name      Value
#
# x = r23   2-over-3 rule
#   = r5sl  exp(slope of the most recent 5 years' log-index)
#
# y = fLBI  Lmean/LF=M
#     fBHE  M/(Z-M) from Beverton-Holt equation
#     fGH   F0.1/(Z-M) from Gedamke-Hoenig
#     fGHe  F0.1/(Z-M) from Gedamke-Hoenig with effort
#
# z = cat3  b = min(1, Icurr/Itrig) for category 3 stocks
#   = cat4  b = 0.8 for category 4 stocks for every fourth year


library(DLMtool)
library(Rcpp)

source('functions/rfb_functions.R')
source('functions/meanlength_fun.R')
sourceCpp('functions/ML_functions.cpp')
source('functions/category3_HCR.R')
source('functions/category4_HCR.R')

source('functions/Pplot2_Blim.R')
source('functions/Tplot_Blim.R')

########### Calculate Blim from MSE
# Blim is the max(0.2 * B0, B such that R/R0 = 0.8)
# Biomass refers to spawning stock biomass
calculate_Blim <- function(steepness, B0, xR = 0.8) {
  h <- steepness
  Blim.SR <- xR * B0 * (1-h) /(4*h - xR * (5*h - 1))
  Blim <- ifelse(Blim.SR < 0.2*B0, 0.2*B0, Blim.SR)
  return(Blim)
}


get_Blim <- function(MSEobj, xR = 0.8, output = c('summary', 'raw')) {
  output <- match.arg(output)
  nm <- MSEobj@nMPs
  nsim <- MSEobj@nsim
  proyears <- MSEobj@proyears
  
  B0 <- MSEobj@OM$SSBMSY/MSEobj@OM$SSBMSY_SSB0
  Blim <- calculate_Blim(steepness = MSEobj@OM$hs, B0 = B0, xR = xR)
  
  if(output == 'summary') {
    PBlim <- matrix(NA, nm, nsim)
    for(m in 1:nm) {
      for(j in 1:nsim) PBlim[m, j] <- sum(MSEobj@SSB[j, m, ] < Blim[j])/proyears * 100
    }
    Blim_BMSY <- mean(Blim/MSEobj@OM$SSBMSY)
    
    MP.summary <- data.frame(MP = MSEobj@MPs, PBlim = round(apply(PBlim, 1, mean, na.rm = T), 2),
                             stdev = round(apply(PBlim, 1, sd, na.rm = T), 2))
    OM.summary <- round(data.frame(Blim_BMSY = mean(Blim/MSEobj@OM$SSBMSY), stdev = sd(Blim/MSEobj@OM$SSBMSY),
                                   Blim_B0 = mean(Blim/B0), stdev = sd(Blim/B0),
                                   BMSY_B0 = mean(MSEobj@OM$SSBMSY_SSB0), stdev = sd(MSEobj@OM$SSBMSY_SSB0)), 2)
    
    return(list(MP.summary = MP.summary, OM.summary = OM.summary))
  }
  if(output == 'raw') {
    B_Blim <- array(NA, dim = c(nsim, nm, proyears))
    for(m in 1:nm) {
      for(j in 1:nsim) B_Blim[j, m, ] <- MSEobj@SSB[j, m, ]/Blim[j]
    }
    return(B_Blim)
  }
}


########### Preliminary MPs
## Update advice every 2 or 3 years
Two_Over_Three_Capped <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind"
  index.samp <- sample_index(x, Data, reps)
  Ind5 <- matrix(index.samp[(nrow(index.samp)-4):nrow(index.samp), ], ncol = reps)
  r <- r_2over3(Ind5, uncertainty_cap = FALSE)
  f <- 1
  b <- 1
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(Two_Over_Three_Capped) <- "Output"
environment(Two_Over_Three_Capped) <- asNamespace("DLMtool")

Islope5y <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind"
  index.samp <- sample_index(x, Data, reps)
  Ind5 <- matrix(index.samp[(nrow(index.samp)-4):nrow(index.samp), ], ncol = reps)
  r <- r_expIslope(Ind5)
  f <- 1
  b <- 1
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(Islope5y) <- "Output"
environment(Islope5y) <- asNamespace("DLMtool")
