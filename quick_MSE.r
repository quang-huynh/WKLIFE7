# Using Version 4.4.1 of DLMtool
# install.packages('DLMtool')

# http://www.datalimitedtoolkit.org/
# https://dlmtool.github.io/DLMtool/userguide/index.html

source('functions/WKLIFE7_HCR.R')

setup()
sfLibrary(DLMtool)
sfExportAll()

OM <- new('OM', Rockfish, Generic_fleet, Precise_Unbiased, Perfect_Imp)

OM@nsim <- 250
MPs <- c('r23_fGH_b3', 'Two_Over_Three_Capped', 'FMSYref', 'NFref')

rockfishMSE <- runMSE(OM=OM, MPs=MPs, interval = 2, reps = 1)

MPs <- c('r23_fLBI_b3', 'r23_fBHE_b3')
OM@nsim <- 250
rockfishMSE2 <- runMSE(OM = OM, MPs = MPs, interval = 2, reps = 1)

save(rockfishMSE, rockfishMSE2, file = 'Rockfish.MSE.RData')
sfStop()

plotFun()

Converge(rockfishMSE)

summary(rockfishMSE)
boxplot(rockfishMSE@B_BMSY[,1,], xlab="Year", ylab="B/BMSY")

Tplot(rockfishMSE)
Tplot2(rockfishMSE)

Pplot(rockfishMSE2)
#Pplot2(rockfishMSE, YVar=c("SSB_SSBMSY", "F_FMSY", "Yield"))
Pplot2(rockfishMSE2, traj = "quant", quants = c(0.2, 0.8), oneIt = FALSE)

VOIplot(rockfishMSE)
VOIplot(rockfishMSE, Par = "OM")
