Tplot_Blim <- function (MSEobj, nam = NA) 
{
  FMSYr <- quantile(MSEobj@F_FMSY, c(0.001, 0.9), na.rm = T)
  BMSYr <- quantile(MSEobj@B_BMSY, c(0.001, 0.975), na.rm = T)
  colsse <- rainbow(100, start = 0, end = 0.36)[1:100]
  colB <- rep(colsse[100], ceiling(BMSYr[2] * 100))
  colB[1:100] <- colsse
  colB <- makeTransparent(colB, 60)
  colsse <- rainbow(200, start = 0, end = 0.36)[200:1]
  colF <- rep(colsse[200], ceiling(FMSYr[2] * 100))
  colF[1:200] <- colsse
  colF <- makeTransparent(colF, 60)
  Yd <- rep(NA, MSEobj@nMPs)
  #PBlim <- rep(NA, MSEobj@nMPS)
  POF <- rep(NA, MSEobj@nMPs)
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  for (mm in 1:MSEobj@nMPs) {
    Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, 
                               mean, na.rm = T)/RefYd, na.rm = T) * 100, 1)
    POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] >= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
                                                                                       mm, ]), na.rm = T) * 100, 1)
  }
  PBlim <- get_Blim(MSEobj, xR = 0.7, output = 'summary')[[1]]$PBlim
  old_par <- par(mfrow = c(2, 1), mai = c(0.9, 1, 0.1, 0.1), 
                 omi = c(0.1, 0.1, 0.4, 0))
  tradeoffplot(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplot(PBlim, Yd, "Prob. biomass < Blim (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 5, hl = 100)
  if (is.na(nam)) 
    mtext(deparse(substitute(MSEobj)), 3, outer = T, line = 0.3, 
          font = 2)
  if (!is.na(nam) & !is.character(nam)) 
    mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
  if (!is.na(nam) & is.character(nam)) 
    mtext(nam, 3, outer = T, line = 0.3, font = 2)
  par(old_par)
}
environment(Tplot_Blim) <- asNamespace('DLMtool')
