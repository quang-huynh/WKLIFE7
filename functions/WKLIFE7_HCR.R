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
#   = r4    1 for category 3 stocks
#
# y = fLBI  Lmean/LF=M
#     fBHE  M/(Z-M) from Beverton-Holt equation
#     fGH   F0.1/(Z-M) from Gedamke-Hoenig
#     fGHe  F0.1/(Z-M) from Gedamke-Hoenig with effort
#
# z = cat3  b = min(1, Icurr/Itrig) for category 3 stocks
#   = cat4  b = 0.8 for category 4 stocks


library(DLMtool)
library(Rcpp)

source('functions/rfb_functions.R')
source('functions/meanlength_fun.R')
sourceCpp('functions/ML_functions.cpp')
source('functions/category3_HCR.R')
source('functions/category4_HCR.R')

########### Preliminary MPs
## Update advice every 2 or 3 years
Two_Over_Three_Capped <- function(x, Data, reps) {
  dependencies = "Data@Cat, Data@CV_Cat, Data@Ind, Data@CV_Ind"
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- r_2over3(Ind5, uncertainty_cap = TRUE)
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
  Ind5 <- sample_index(x, Data, reps, nyrs = 5)
  r <- r_expIslope(Ind5)
  f <- 1
  b <- 1
  Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
  TAC <- Cc * r * f * b
  TACfilter(TAC)
}
class(Islope5y) <- "Output"
environment(Islope5y) <- asNamespace("DLMtool")
