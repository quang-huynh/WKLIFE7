Pplot2_Blim <- function (MSEobj, Hist, YVar = c("SSB_SSBlim", "F_FMSY", "Yield"), MPs = NA, 
          sims = NULL, traj = c("all", "quant"), quants = c(0.1, 0.9), 
          incquant = TRUE, quantcol = "lightgray", RefYield = c("lto", "curr"), 
          HistoricalYr = FALSE, LastYr = TRUE, maxMP = 6, alpha = 60, cex.axis = 1.35, 
          cex.lab = 1.4, YLab = NULL, incMP = TRUE, MPcex = 1.4, incLeg = TRUE, 
          cex.leg = 1.5, legPos = "topleft", yline = NULL, parOR = FALSE, 
          xaxis = TRUE, yaxis = TRUE, oneIt = TRUE, ...) 
{
  YVars <- c("SSB_SSB0", "SSB_SSBMSY", "F_FMSY", "Yield", "SSB_SSBlim")  # Add for Blim
  YVar <- match.arg(YVar, choices = YVars, several.ok = TRUE)
  op <- par(no.readonly = TRUE)
  if (!is.null(YLab) & length(YLab) != length(YVar)) 
    stop("Length of YLab must equal length of YVar")
  if (!is.null(sims) & all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims)
  if (!is.null(sims) & all(!is.na(MPs))) 
    MSEobj <- Sub(MSEobj, sims = sims, MPs = MPs)
  if (is.null(sims) & !all(is.na(MPs))) 
    MSEobj <- Sub(MSEobj, MPs = MPs)
  nsim <- MSEobj@nsim
  if (nsim < 10) 
    alpha <- 180
  nMPs <- MSEobj@nMPs
  if (nMPs > maxMP) {
    message("MSE object has more than ", maxMP, " MPs. Plotting the first ", 
            maxMP)
    MSEobj <- Sub(MSEobj, MPs = 1:maxMP)
    nMPs <- MSEobj@nMPs
  }
  MPs <- MSEobj@MPs
  proyears <- MSEobj@proyears
  nyears <- MSEobj@nyears
  RefYd <- MSEobj@OM$RefY
  RefYield <- match.arg(RefYield)
  traj <- match.arg(traj)
  temp <- as.matrix(expand.grid(1:nsim, 1:nMPs, 1:proyears))
  Deplet <- array(NA, dim = dim(MSEobj@B_BMSY))
  Deplet[temp] <- (MSEobj@B_BMSY[temp] * MSEobj@OM$SSBMSY_SSB0[temp[, 
                                                                    1]])
  pastC <- apply(MSEobj@CB_hist[, , , , drop = FALSE], c(1, 
                                                         3), sum, na.rm = TRUE)/RefYd
  temp <- aperm(replicate(nMPs, pastC), c(1, 3, 2))
  lastYr <- temp[, , MSEobj@nyears, drop = FALSE]
  Yield <- abind::abind(lastYr, MSEobj@C[, , , drop = FALSE]/RefYd, 
                        along = 3)
  SSB_SSBlim <- get_Blim(MSEobj, xR = 0.7, output = 'raw') # Summary Blim
  Dat <- list(SSB_SSB0 = Deplet, SSB_SSBMSY = MSEobj@B_BMSY, 
              F_FMSY = MSEobj@F_FMSY, Yield = Yield, SSB_SSBlim = SSB_SSBlim)
  Dat <- Dat[YVar]
  if ("Yield" %in% YVar & RefYield == "curr") {
    ny <- dim(Dat$Yield)[3]
    Dat$Yield <- Dat$Yield[, , , drop = FALSE]/Dat$Yield[, 
                                                         , rep(1, ny), drop = FALSE]
  }
  if ("Yield" %in% YVar & !LastYr) {
    Dat$Yield <- Dat$Yield[, , 2:proyears, drop = FALSE]
  }
  
  # Historical years
  SSB_hist <- apply(MSEobj@SSB_hist, c(1, 3), sum)
  FM_hist <- Hist$SampPars$qs[sims] * Hist$SampPars$Find[sim, ]
  CB_hist <- apply(MSEobj@CB_hist, c(1, 3), sum) # sum over areas and ages
  if(RefYield == "lto") Yield.hist <- CB_hist/MSEobj@OM$RefY
  if(RefYield == "curr") Yield.hist <- CB_hist/CB_hist[, nyears]
  
  Blim.absolute <- get_Blim(MSEobj, xR = 0.7, output = 'raw', magnitude = 'absolute')
  DatHist <- list(SSB_SSB0 = SSB_hist/MSEobj@OM$SSB0,
                  SSB_SSBMSY = SSB_hist/MSEobj@OM$SSBMSY,
                  F_FMSY = FM_hist/MSEobj@OM$FMSY_M,
                  Yield = Yield.hist,
                  SSB_SSBlim = SSB_hist/Blim.absolute)
  DatHist <- DatHist[YVar]
  
  nr <- length(Dat)
  nc <- nMPs
  dots <- list(...)
  ylims <- cbind(0, unlist(lapply(Dat, quantile, 0.9, na.rm = TRUE)))
  if ("SSB_SSB0" %in% YVar) {
    index <- which(YVar == "SSB_SSB0")
    ylims[index, ] <- c(0, max(1, max(ylims[index, ])))
  }
  if (length(dots$ylim) != 0) 
    ylims <- matrix(rep((dots$ylim), length(Dat)), nrow = nr, 
                    byrow = TRUE)
  colrange <- matrix(unlist(lapply(Dat, quantile, c(0.001, 
                                                    0.975), na.rm = TRUE)), nrow = nr, byrow = TRUE)
  colsse <- rainbow(100, start = 0, end = 0.36)[1:100]
  Col <- makeTransparent(colsse, alpha)
  if (length(dots$lwd) == 0) 
    lwd <- 3
  if (length(dots$lwd) != 0) 
    lwd <- dots$lwd
  YLabs <- list(expression(SSB/SSB[0]), expression(SSB/SSB[MSY]), 
                expression(F/F[MSY]), "Yield relative\n to Long-Term\n Optimum", expression(SSB/SSB[lim]))
  if ("Yield" %in% YVar & RefYield == "curr") 
    YLabs[[4]] <- expression(Yield/Yield[current])
  YLabs <- YLabs[match(YVar, YVars)]
  if (!is.null(YLab)) 
    YLabs <- YLab
  if (!parOR) {
    if ("Yield" %in% YVar & RefYield != "curr") {
      op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 
                                                      2, 0, 0), oma = c(4, 8, 2, 1))
    }
    else op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2, 
                                                         2, 0, 0), oma = c(4, 4, 2, 1))
  }
  if (parOR) {
    nr <- par()$mfrow[1]
    nc <- par()$mfrow[2]
  }
  for (X in 1:length(Dat)) {
    Col2 <- Col
    dat <- Dat[[X]]
    datHist <- DatHist[[X]]
    ylim <- ylims[X, ]
    ylab <- YLabs[[X]]
    if (grepl("F_FMSY", YVar[X])) 
      Col2 <- rev(Col)
    for (mm in 1:nMPs) {
      if(!HistoricalYr) x.lower <- 1
      if(HistoricalYr) x.lower <- -nyears
      x.upper <- length(dat[1, mm, ])
      plot(1:length(dat[1, mm, ]), dat[1, mm, ], ylim = ylim, xlim = c(x.lower, x.upper),
           type = "n", axes = FALSE, xlab = "", ylab = "")
      if(HistoricalYr) abline(v = 0, lwd = 1, lty = 3)
      if (traj == "all")  
        for (i in 1:MSEobj@nsim) {
          lines(dat[i, mm, ], col = Col2[min(100, ceiling(dat[i, mm, length(dat[i, mm, ])] * 100))], lwd = lwd)
          if(HistoricalYr) lines(x.lower:-1, datHist[i, ], 
                                 col = Col2[min(100, ceiling(dat[i, mm, length(dat[i, mm, ])] * 100))], lwd = lwd)
        }
      if (traj == "quant") {
        stats <- apply(dat[, mm, , drop = FALSE], 3, 
                       quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
        statsHist <- apply(datHist, 2, quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
        if (length(quants) == 4) {
          stats2 <- apply(dat[, mm, , drop = FALSE], 
                          3, quantile, c(quants[3], quants[4]), na.rm = TRUE)
          stats2Hist <- apply(datHist, 3, quantile, c(quants[3], quants[4]), na.rm = TRUE)
        }
        if (length(quants) != 4) 
          stats2 <- stats2Hist <- NULL
        if (!incquant) {
          lines(1:length(stats[2, ]), stats[2, ], lwd = 3)
          if(HistoricalYr) lines(x.lower:-1, statsHist[2, ], lwd = 3)
        }
        if (incquant) {
          if (!is.null(stats2)) {
            if (length(quantcol) < 2) {
              quantcol <- c(quantcol, "darkgray")
              polygon(x = c(1:length(stats2[2, ]), length(stats2[2, ]):1), 
                      y = c(stats2[1, ], rev(stats2[2, ])), col = quantcol[2], border = FALSE)
            }
          }
          polygon(x = c(1:length(stats[2, ]), length(stats[2, ]):1), y = c(stats[1, ], rev(stats[3, ])), 
                  col = quantcol[1], border = FALSE)          
          lines(1:length(stats[2, ]), stats[2, ], lwd = 3)
          if(HistoricalYr) {
            polygon(x = c(-nyears:-1, -1:-nyears), y = c(statsHist[1, ], rev(statsHist[3, ])),
                    col = quantcol[1], border = FALSE)
            lines(-nyears:-1, statsHist[2, ], lwd = 3)
            if(oneIt) 
              lines(-nyears:-1, datHist[2, ], lwd = 1)
          }
          if (oneIt) 
            lines(dat[2, mm, , drop = FALSE], lwd = 1)
        }
      }
      if (X == nr & xaxis) 
        axis(side = 1, labels = TRUE, cex.axis = cex.axis)
      if (X != nr | !xaxis) 
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
      if (incMP & X == 1 & !parOR) 
        mtext(side = 3, MSEobj@MPs[mm], cex = MPcex)
      if (incMP & parOR) 
        mtext(side = 3, MSEobj@MPs[mm], cex = MPcex)
      if (mm == 1 & incLeg & traj == "quant" & X == 1) {
        if (incquant) {
          if (is.null(stats2)) {
            pquants <- quants
            if (max(quants) < 1) 
              pquants <- quants * 100
            legend(legPos, legend = c(expression("50"^"th" * 
                                                   " (median)"), bquote(.(pquants[1])^th ~ 
                                                                          and ~ .(pquants[2])^th)), pch = 22, title = "Percentile", 
                   pt.bg = c("black", quantcol[1], quantcol[2]), 
                   bty = "n", cex = cex.leg, xpd = NA, col = "black")
          }
          if (!is.null(stats2)) {
            pquants <- quants
            if (max(quants) < 1) 
              pquants <- quants * 100
            legend(legPos, legend = c(expression("50"^"th" * 
                                                   " (median)"), bquote(.(pquants[1])^th ~ 
                                                                          and ~ .(pquants[2])^th), bquote(.(pquants[3])^th ~ 
                                                                                                            and ~ .(pquants[4])^th)), pch = 22, title = "Percentile", 
                   pt.bg = c("black", quantcol[1], quantcol[2]), 
                   bty = "n", cex = cex.leg, xpd = NA, col = "black")
          }
        }
        else {
          legend(legPos, legend = "Median", pch = 15, 
                 col = "black", bty = "n", cex = 1.25)
        }
      }
      if (!is.null(yline)) 
        abline(h = yline[X], lwd = 2, col = "darkgray", 
               lty = 2)
    }
  }
  mtext(side = 1, "Projection Years", line = 2, cex = cex.lab, 
        outer = TRUE)
  if (!parOR) 
    par(op)
  invisible(Dat)
}
