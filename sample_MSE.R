# Using Version 4.4.1 of DLMtool
# install.packages('DLMtool')

# http://www.datalimitedtoolkit.org/
# https://dlmtool.github.io/DLMtool/userguide/index.html

# Load OMs
load('basicOM.RData')

# Load HCR
source('functions/WKLIFE7_HCR.R')

category3.HCR <- paste0(rep(c('r23_', 'r5sl_'), each = 4), 
                        c('fLBI_', 'fBHE_', 'fGH_', 'fGHe_'), 'cat3')

category4.HCR <- paste0(c('fLBI_', 'fBHE_', 'fGH_'), 'cat4')


setup()
sfLibrary(DLMtool)
sfExportAll()
clusterEvalQ(sfGetCluster(), expr = {
  Rcpp::sourceCpp(file = "functions/ML_functions.cpp")
}
)


### RUN MSE
haddockMSE <- runMSE(OM = haddockOM, MPs = c('FMSYref', 'NFref', category3.HCR), 
                     interval = 2, reps = 1)
#save(haddockMSE, file = "haddockMSE.RData")
sfStop()

## Summary statistics
summary(haddockMSE)
get_Blim(haddockMSE, xR = 0.7) # R/R0 = 0.7

# Plot MSE figures

Pplot2_Blim(haddockMSE, YVar = c('SSB_SSBlim', 'F_FMSY'), maxMP = 10,
            traj = "quant", quants = c(0.1, 0.9), oneIt = FALSE)

plotFun()

Tplot(haddockMSE)

Pplot(haddockMSE)
VOIplot(haddockMSE)
VOIplot(haddockMSE, Par = "OM")

VOI(haddockMSE)
VOI2(haddockMSE)
