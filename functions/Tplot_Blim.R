Tplot_Blim <- function (MSEobj, nam = NA) 
{
  
  tradeoffplotQCH <- function (x, y, xlab, ylab, labs, cex, vl, hl) # Stolen  from Tobias
  {
    adjj <- c(0.7, 1.3)
    XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, 
                                                               na.rm = T) * adjj, 110)))
    YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, 
                                                               na.rm = T) * adjj, 110)))
    coly <- gplots::rich.colors(length(x))
    #coly <- rep(brewer.pal(8,"Dark2"), 20)[1:length(labs)]
    #coly[labs %in% c("AvC", "curE", "FMSYref")] <- "black"
    plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
    abline(v = vl, col = "#99999940", lwd = 2)
    abline(h = hl, col = "#99999940", lwd = 2)
    points(x, y, col = coly, pch=16, cex=1.2)
    text(x, y, labels = 1:length(x), pos = 3)
    segments(x0=x,x1=x,y0=-30,y1=y, col=coly, lty=2)
    segments(x0=-30,x1=x,y0=y,y1=y, col=coly, lty=2)    
  }
  
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
  P10 <- rep(NA, MSEobj@nMPs)
  P50 <- rep(NA, MSEobj@nMPs)
  #PBlim <- rep(NA, MSEobj@nMPS)
  POF <- rep(NA, MSEobj@nMPs)
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  RefYd <- MSEobj@OM$RefY
  for (mm in 1:MSEobj@nMPs) {
    Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, 
                               mean, na.rm = T)/RefYd, na.rm = T) * 100, 1)
    POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] >= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
                                                                                       mm, ]), na.rm = T) * 100, 1)
    P10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
    P50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
                                                                                         mm, ])) * 100, 1)
  }
  PBlim <- get_Blim(MSEobj, xR = 0.7, output = 'summary')[[1]]$PBlim
  browser()
  old_par <- par(mfrow = c(2, 1), mar = c(5, 5, 1, 1), 
                 oma = c(0, 0, 4, 0))
  layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,5,5,5,5), ncol = 5))
  tradeoffplotQCH(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplotQCH(PBlim, Yd, "Prob. biomass < Blim (%)", "Relative yield", 
               MSEobj@MPs[1:MSEobj@nMPs], vl = 5, hl = 100)
  tradeoffplotQCH(P50, Yd, "Prob. biomass < 0.5BMSY (%)", "Relative yield", 
                  MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  tradeoffplotQCH(P10, Yd, "Prob. biomass < 0.1BMSY (%)", "Relative yield", 
                  MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
  
  par(mar = c(5, 0.5, 1, 0.5))
  plot(1, type="n", axes=FALSE, xlab="", ylab="")    
  legend("left", legend = paste0(MSEobj@MPs[1:MSEobj@nMPs], " (", 1:MSEobj@nMPs, ")"),
         col = gplots::rich.colors(length(MSEobj@MPs)), bty="n",pch=16, cex=0.9)    
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
