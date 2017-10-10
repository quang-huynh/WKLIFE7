
get_Pplot_median <- function(MSEobj, Hist, RefYield = c('lto', 'curr'), LastYr = F, MPs = NA,
                             sims = NULL) {
  if (!is.null(sims) & all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims)
  if (!is.null(sims) & all(!is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
  if (is.null(sims) & !all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, MPs = MPs)
  nsim <- MSEobj@nsim
  if (nsim < 10) alpha <- 180
  nMPs <- MSEobj@nMPs
  MPs <- MSEobj@MPs
  proyears <- MSEobj@proyears
  nyears <- MSEobj@nyears
  RefYd <- MSEobj@OM$RefY
  RefYield <- match.arg(RefYield)
  temp <- as.matrix(expand.grid(1:nsim, 1:nMPs, 1:proyears))
  Deplet <- array(NA, dim = dim(MSEobj@B_BMSY))
  Deplet[temp] <- (MSEobj@B_BMSY[temp] * MSEobj@OM$SSBMSY_SSB0[temp[, 1]])
  pastC <- apply(MSEobj@CB_hist[, , , , drop = FALSE], c(1, 
                                                         3), sum, na.rm = TRUE)/RefYd
  temp <- aperm(replicate(nMPs, pastC), c(1, 3, 2))
  lastYr <- temp[, , MSEobj@nyears, drop = FALSE]
  Yield <- abind::abind(lastYr, MSEobj@C[, , , drop = FALSE]/RefYd, 
                        along = 3)
  SSB_SSBlim <- get_Blim(MSEobj, xR = 0.7, output = 'raw') # Summary Blim
  Dat <- list(SSB_SSB0 = Deplet, SSB_SSBMSY = MSEobj@B_BMSY, 
              F_FMSY = MSEobj@F_FMSY, Yield = Yield, SSB_SSBlim = SSB_SSBlim)
  if (RefYield == "curr") {
    ny <- dim(Dat$Yield)[3]
    Dat$Yield <- Dat$Yield[, , , drop = FALSE]/Dat$Yield[, 
                                                         , rep(1, ny), drop = FALSE]
  }
  if (!LastYr) {
    Dat$Yield <- Dat$Yield[, , 2:proyears, drop = FALSE]
  }
  
  # Historical years
  SSB_hist <- apply(MSEobj@SSB_hist, c(1, 3), sum)
  FM_hist <- Hist$SampPars$qs[1:nsim] * Hist$SampPars$Find[1:nsim, ]
  CB_hist <- apply(MSEobj@CB_hist, c(1, 3), sum) # sum over areas and ages
  if(RefYield == "lto") Yield.hist <- CB_hist/MSEobj@OM$RefY
  if(RefYield == "curr") Yield.hist <- CB_hist/CB_hist[, nyears]
  
  
  Blim.absolute <- get_Blim(MSEobj, xR = 0.7, output = 'raw', magnitude = 'absolute')
  DatHist <- list(SSB_SSB0 = SSB_hist/MSEobj@OM$SSB0,
                  SSB_SSBMSY = SSB_hist/MSEobj@OM$SSBMSY,
                  F_FMSY = FM_hist/MSEobj@OM$FMSY,
                  Yield = Yield.hist,
                  SSB_SSBlim = SSB_hist/Blim.absolute)
  
  Datmed <- lapply(Dat, function(x) apply(x, c(2, 3), median, na.rm = TRUE))
  DatHistmed <- lapply(DatHist, function(x) apply(x, 2, median, na.rm = TRUE))
  
  output <- list(Dat = Datmed, DatHist = DatHistmed, MPs = MPs)
  return(output)
}


Pplot2_Blim_median <- function(
  MSE.list, YVar = c("SSB_SSB0", "F_FMSY", "Yield"), Stocks = Stock.names[ind],
  RefYield = c('lto', 'curr'), MPs = MSE.list[[1]]$MPs, 
  col.MSE = gplots::rich.colors(length(MSE.list)),
  HistoricalYr = FALSE, LastYr = TRUE, alpha = 60, cex.axis = 1.35, 
  cex.lab = 1.4, YLab = NULL, incMP = TRUE, MPcex = 1.4, cex.stock = 1.4,
  cex.leg = 1.5, legPos = "topleft", yline = NULL, parOR = FALSE, 
  xaxis = TRUE, yaxis = TRUE, ExtraYrs = 5, legendWidth = 1, ...) 
{
  YVars <- c("SSB_SSB0", "SSB_SSBMSY", "F_FMSY", "Yield", "SSB_SSBlim")  # Add for Blim
  YVar <- match.arg(YVar, choices = YVars, several.ok = TRUE)
  RefYield <- match.arg(RefYield)
  op <- par(no.readonly = TRUE)
  nr <- length(YVar)
  nc <- length(MSE.list[[1]]$MPs)
  dots <- list(...)
  
  Dat1 <- do.call(cbind, lapply(MSE.list, function(x) x$Dat[[1]]))
  Dat2 <- do.call(cbind, lapply(MSE.list, function(x) x$Dat[[2]]))
  Dat3 <- do.call(cbind, lapply(MSE.list, function(x) x$Dat[[3]]))
  Dat4 <- do.call(cbind, lapply(MSE.list, function(x) x$Dat[[4]]))
  Dat5 <- do.call(cbind, lapply(MSE.list, function(x) x$Dat[[5]]))
  Dat.all <- list(SSB_SSB0 = Dat1, SSB_SSBMSY = Dat2, F_FMSY = Dat3, 
                  Yield = Dat4, SSB_SSBlim = Dat5)
  
  ylims <- cbind(0, unlist(lapply(Dat.all, quantile, 0.9, na.rm = TRUE)))
  if ("SSB_SSB0" %in% YVar) {
    index <- which(YVar == "SSB_SSB0")
    ylims[index, ] <- c(0, max(1, max(ylims[index, ])))
  }
  if (length(dots$ylim) != 0) 
    ylims <- matrix(rep((dots$ylim), length(Dat.all)), nrow = nr, 
                    byrow = TRUE)
  
  if (length(dots$lwd) == 0) 
    lwd <- 3
  if (length(dots$lwd) != 0) 
    lwd <- dots$lwd
  YLabs <- list(expression(SSB/SSB[0]), expression(SSB/SSB[MSY]), 
                expression(F/F[MSY]), "Yield relative\n to Long-Term\n Optimum", expression(SSB/SSB[lim]))
  if (RefYield == "curr") 
    YLabs[[4]] <- expression(Yield/Yield[current])
  if (!is.null(YLab)) 
    YLabs <- YLab
  if (!parOR) {
    if ("Yield" %in% YVar & RefYield != "curr") {
      op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 2, 0, 0), oma = c(4, 8, 2, 1))
    }
    else op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 2, 0, 0), 
                   oma = c(4, 4, 2, 1))
    layout(cbind(matrix(1:(nr * nc), ncol = nc, byrow = TRUE), 
                matrix(nr*nc+1, ncol = 1, nrow = nr)), widths = c(rep(1, nc), legendWidth))
  }
  if (parOR) {
    nr <- par()$mfrow[1]
    nc <- par()$mfrow[2]
    layout(cbind(matrix(1:(nr * nc), ncol = nc, byrow = TRUE), 
                 matrix(nr*nc+1, ncol = 1, nrow = nr)), widths = c(rep(1, nc), legendWidth))
  }
  loop <- match(YVar, YVars)
  for (X in loop) { # Loop over data-type
    #ylim <- ylims[X, ]
    if(YVars[X] == 'SSB_SSB0') ylim <- c(0, 1.1)
    if(YVars[X] == 'SSB_SSBMSY') ylim <- c(0, 5)
    if(YVars[X] == 'F_FMSY') ylim <- c(0, 1.5)
    if(YVars[X] == 'Yield') ylim <- c(0, 1.5)
    if(YVars[X] == 'SSB_SSBlim') ylim <- c(0, 4)
    
    ylab <- YLabs[[X]]
    for (mm in 1:length(MPs)) { # Loop over MPs
      if(!HistoricalYr) x.lower <- 1
      if(HistoricalYr) x.lower <- -length(MSE.list[[1]]$DatHist[[1]])
      x.upper <- dim(MSE.list[[1]]$Dat[[1]])[2]
      
      # First OM for dimensions
      dat <- MSE.list[[1]]$Dat[[X]] # dat is a matrix of mp rows and year cols
      datHist <- MSE.list[[1]]$DatHist[[X]]
      plot(1:length(dat[mm, ]), dat[mm, ], ylim = ylim, xlim = c(x.lower, x.upper + ExtraYrs),
           type = "n", axes = FALSE, xlab = "", ylab = "")
      if(HistoricalYr) abline(v = 0, lwd = 1, lty = 3)
      
      for (i in 1:length(MSE.list)) { # Loop over all OMs
        dat <- MSE.list[[i]]$Dat[[X]]
        datHist <- MSE.list[[i]]$DatHist[[X]]
        lines(dat[mm, ], col = col.MSE[i], lwd = lwd)
        text(length(dat[mm, ]), dat[mm, length(dat[mm, ])], labels = paste0('(', i, ')'), pos = 4,
             cex = cex.stock)
        if(HistoricalYr) lines(x.lower:-1, datHist, col = col.MSE[i], lwd = lwd)
        if(YVars[X] != 'SSB_SSB0') abline(h = 1, lty = 2)
      }
      
      if (X == max(match(YVar, YVars)) & xaxis) 
        axis(side = 1, labels = TRUE, cex.axis = cex.axis)
      if (X != max(match(YVar, YVars)) | !xaxis) 
        axis(side = 1, labels = FALSE)
      if (mm == 1 & yaxis) {
        axis(side = 2, labels = TRUE, cex.axis = cex.axis, 
             las = 1)
        mtext(side = 2, ylab, cex = cex.lab, line = 3)
      }
      if (parOR) {
        if (xaxis) 
          axis(side = 1, labels = TRUE, cex.axis = cex.axis, 
               las = 1)
        if (yaxis & mm == 1) 
          axis(side = 2, labels = TRUE, cex.axis = cex.axis, 
               las = 1)
      }
      axis(side = 2, labels = FALSE)
      if (incMP & X == loop[1] & !parOR) 
        mtext(side = 3, MPs[mm], cex = MPcex)
      if (incMP & parOR) 
        mtext(side = 3, MPs[mm], cex = MPcex)
      
      if (!is.null(yline)) 
        abline(h = yline[X], lwd = 2, col = "darkgray", 
               lty = 2)
    }
  }
  # Legend plot
  plot(1:10, 1:10, axes = FALSE, ylab = '', xlab = '', type = 'n')
  legend(legPos, paste0(Stocks, ' (', c(1:length(Stocks)), ')'), col = col.MSE, 
         pch = 16, lty = 1, lwd = 3, bty = 'n', cex = cex.stock)
  
  mtext(side = 1, "Projection Years", line = 2, cex = cex.lab, 
        outer = TRUE)
  if (!parOR) 
    par(op)
  invisible()
}
