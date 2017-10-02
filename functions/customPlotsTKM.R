

require(RColorBrewer)

### not happy with some default settings in the DLMtool plots + additional graphs



tradeoffplotTKM <- function (x, y, xlab, ylab, labs, cex, vl, hl) 
{
    adjj <- c(0.7, 1.3)
    XLim <- c(min(c(-10, min(x, na.rm = T) * adjj)), max(c(max(x, 
        na.rm = T) * adjj, 110)))
    YLim <- c(min(c(-10, min(y, na.rm = T) * adjj)), max(c(max(y, 
        na.rm = T) * adjj, 110)))
    coly <- rep(brewer.pal(8,"Dark2"), 20)[1:length(labs)]
    coly[labs %in% c("AvC", "curE", "FMSYref")] <- "black"
    plot(NA, xlim = XLim, ylim = YLim, xlab = xlab, ylab = ylab)
    abline(v = vl, col = "#99999940", lwd = 2)
    abline(h = hl, col = "#99999940", lwd = 2)
    points(x, y, col = coly, pch=16, cex=1.2)
    segments(x0=x,x1=x,y0=-30,y1=y, col=coly, lty=2)
    segments(x0=-30,x1=x,y0=y,y1=y, col=coly, lty=2)    
}



TplotTKM <- function (MSEobj, nam = NA, main = NA) 
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
    P10 <- rep(NA, MSEobj@nMPs)
    P50 <- rep(NA, MSEobj@nMPs)
    P100 <- rep(NA, MSEobj@nMPs)
    POF <- rep(NA, MSEobj@nMPs)
    yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
    RefYd <- MSEobj@OM$RefY
    for (mm in 1:MSEobj@nMPs) {
        Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, mean, 
            na.rm = T)/RefYd, na.rm = T) * 100, 1)
        POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] >= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
            mm, ]), na.rm = T) * 100, 1)
        P10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
        P50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
        P100[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
    }
    coly <- rep(brewer.pal(8,"Dark2"), 20)[1:length(MSEobj@MPs[1:MSEobj@nMPs])]    
    coly[MSEobj@MPs[1:MSEobj@nMPs] %in% c("AvC", "curE", "FMSYref")] <- "black"    
    old_par <- par(mfrow = c(2, 3), mai = c(0.9, 0.8, 0.1, 0.1), 
        omi = c(0.1, 0.1, 0.4, 0))
    tradeoffplotTKM(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
                 MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
    tradeoffplotTKM(P100, Yd, "Prob. biomass < BMSY (%)", "Relative yield", 
                    MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")    
    tradeoffplotTKM(P50, Yd, "Prob. biomass < 0.5BMSY (%)", "Relative yield", 
                    MSEobj@MPs[1:MSEobj@nMPs], vl = c(5,50), hl = 100)
    tradeoffplotTKM(P10, Yd, "Prob. biomass < 0.1BMSY (%)", "Relative yield", 
                    MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")    
    legend("bottomleft", legend=MSEobj@MPs[1:MSEobj@nMPs],
           col = coly, bty="n",pch=16, cex=0.9)    
    if (is.na(nam)) 
        ##        mtext(deparse(substitute(MSEobj)), 3, outer = T, line = 0.3, font = 2)
        mtext(main, 3, outer = T, line = 0.3, font = 2)
    if (!is.na(nam) & !is.character(nam)) 
        mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
    if (!is.na(nam) & is.character(nam)) 
        mtext(nam, 3, outer = T, line = 0.3, font = 2)
    par(old_par)
}



Tplot2TKM <- function (MSEobj, nam = NA, main=NA) 
{
    LTY <- rep(NA, MSEobj@nMPs)
    STY <- rep(NA, MSEobj@nMPs)
    VY <- rep(NA, MSEobj@nMPs)
    B10 <- rep(NA, MSEobj@nMPs)
    yend <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
    ystart <- 1:5
    RefYd <- MSEobj@OM$RefY
    y1 <- 1:(MSEobj@proyears - 1)
    y2 <- 2:MSEobj@proyears
    for (mm in 1:MSEobj@nMPs) {
        LTY[mm] <- round(sum(MSEobj@C[, mm, yend]/RefYd > 0.5, 
            na.rm = T)/(MSEobj@nsim * length(yend)), 3) * 100
        STY[mm] <- round(sum(MSEobj@C[, mm, ystart]/RefYd > 0.5, 
            na.rm = T)/(MSEobj@nsim * length(ystart)), 3) * 100
        AAVY <- apply(((MSEobj@C[, mm, y1] - MSEobj@C[, mm, y2])^2)^0.5, 
            1, mean, na.rm = T)/apply(MSEobj@C[, mm, y2], 1, 
            mean, na.rm = T)
        VY[mm] <- round(sum(AAVY < 0.1, na.rm = T)/MSEobj@nsim, 
            3) * 100
        B10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] > 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])), 3) * 100
    }
    op <- par(mfrow = c(1, 2), mai = c(1.5, 1.5, 0.1, 0.1), omi = c(0.1, 
        0.1, 0.4, 0.4))
    tradeoffplotTKM(STY, LTY, "P(Short term yield over half FMSY)", 
        "P(Long term yield over half FMSY)", MSEobj@MPs[1:MSEobj@nMPs], 
        vl = 1, hl = 1)
    tradeoffplotTKM(B10, VY, "P(Biomass > 0.1 BMSY)", "P(CV in yield less than 0.1)", 
        MSEobj@MPs[1:MSEobj@nMPs], vl = 1, hl = 1)
    if (is.na(nam)) 
        mtext(main, 3, outer = T, line = 5, font = 2)   
    if (!is.na(nam)) 
        mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
    par(op)
}



HplotTKM <- function(MSEobj, main=NA,
                     quants = c(0.1,0.9),
                     parOR = FALSE){

    ny <- MSEobj@nyears
    py <- MSEobj@proyears
    nc <- MSEobj@nMPs

    ## SSB
    ssbhist <- apply(resMSE@SSB_hist, c(1,3), sum, na.rm=TRUE) ##sum over ages and area
    statsHist <- apply(ssbhist, 2,
                           quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
    statsSSBmsy <- quantile(MSEobj@OM$SSBMSY, probs= c(quants[1], 0.5, quants[2]), na.rm = TRUE)

    ## Fishing mortality
    statsHistF <- apply(resMSE@FM_hist, 3,
                           quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
    statsFFmsy <- quantile(MSEobj@OM$FMSY, probs= c(quants[1], 0.5, quants[2]), na.rm = TRUE)
    
    ylim <- range(MSEobj@SSB, MSEobj@SSB_hist)
    ylimF <- range(MSEobj@FM, MSEobj@FM_hist)    

    nr <- 2  ## 2 if adding F
    if (!parOR) {
        op <- par(mfrow = c(nr, nc), bty = "n", mar = c(2,2, 0, 0), oma = c(4, 4, 8, 1))
    }
    if (parOR) {
        nr <- par()$mfrow[1]
        nc <- par()$mfrow[2]
    }
          
    for (mm in 1:nc){
        statsPro <- apply(MSEobj@SSB[, mm, , drop = FALSE], 3,
                      quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
        
        plot(1:(ny+py), rep(NA,(ny+py)),
             ylim=ylim, xlab="", ylab="")
        mtext(MSEobj@MPs[mm],3,1, font=2)
        if(mm==1) mtext("SSB",2,3)        

        polygon(x=c(1:(ny+py), (ny+py):1),
                y=c(rep(statsSSBmsy[1],(ny+py)),rep(statsSSBmsy[3],(ny+py))),
                col=rgb(192/255,192/255,192/255,0.9), border=FALSE)
        lines(1:(ny+py), rep(statsSSBmsy[2],(ny+py)), col="grey40",lty=2, lwd=2)
        
        abline(v=ny+0.5, lty=2, col="darkorange")
        
        polygon(x = c(1:length(statsHist[2, ]), length(statsHist[2,]):1),
                y = c(statsHist[1, ], rev(statsHist[3, ])),
                col=rgb(0/255,102/255,204/255,0.4), border = FALSE)
        polygon(x = c((ny+1):(ny+length(statsPro[2, ])), (ny+length(statsPro[2,])):(ny+1)),
                y = c(statsPro[1, ], rev(statsPro[3, ])),
                col=rgb(0/255,102/255,204/255,0.4), border = FALSE)
        
        lines(1:length(statsHist[2, ]), statsHist[2, ], lwd = 2,
              col=rgb(0/255,51/255,102/255,0.9))              
        lines((ny+1):(length(statsPro[2, ])+ny), statsPro[2, ], lwd = 2,
              col=rgb(0/255,51/255,102/255,0.9))
    }

            
    for (mm in 1:nc){
        statsProF <- apply(MSEobj@FM[, mm, , drop = FALSE], 3,
                      quantile, c(quants[1], 0.5, quants[2]), na.rm = TRUE)
        
        plot(1:(ny+py), rep(NA,(ny+py)),
             ylim=ylimF, xlab="", ylab="")
        mtext("Years",1,3)
        if(mm==1) mtext("Fishing mortality",2,3)
    
        polygon(x=c(1:(ny+py), (ny+py):1),
                y=c(rep(statsFFmsy[1],(ny+py)),rep(statsFFmsy[3],(ny+py))),
                col=rgb(192/255,192/255,192/255,0.9), border=FALSE)
        lines(1:(ny+py), rep(statsFFmsy[2],(ny+py)), col="grey40",lty=2, lwd=2)
        
        abline(v=ny+0.5, lty=2, col="darkorange")
        
        polygon(x = c(1:length(statsHistF[2, ]), length(statsHistF[2,]):1),
                y = c(statsHistF[1, ], rev(statsHistF[3, ])),
                col=rgb(0/255,153/255,0/255,0.4), border = FALSE)
        polygon(x = c((ny+1):(ny+length(statsProF[2, ])), (ny+length(statsProF[2,])):(ny+1)),
                y = c(statsProF[1, ], rev(statsProF[3, ])),
                col=rgb(0/255,153/255,0/255,0.4), border = FALSE)
        
        lines(1:length(statsHistF[2, ]), statsHistF[2, ], lwd = 2,
              col=rgb(0/255,51/255,0/255,0.9))              
        lines((ny+1):(length(statsProF[2, ])+ny), statsProF[2, ], lwd = 2,
              col=rgb(0/255,51/255,0/255,0.9))
    }
        mtext(main, 3, outer = T, line = 5, font = 2)    
}





TplotsingleTKM <- function (MSEobj, nam = NA, main = NA) 
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
    P10 <- rep(NA, MSEobj@nMPs)
    P50 <- rep(NA, MSEobj@nMPs)
    P100 <- rep(NA, MSEobj@nMPs)
    POF <- rep(NA, MSEobj@nMPs)
    yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
    RefYd <- MSEobj@OM$RefY
    for (mm in 1:5){     ##1:MSEobj@nMPs) {
        
        Yd[mm] <- round(mean(apply(MSEobj@C[, mm, yind], 1, mean, 
            na.rm = T)/RefYd, na.rm = T) * 100, 1)
        POF[mm] <- round(sum(MSEobj@F_FMSY[, mm, ] >= 1, na.rm = T)/prod(dim(MSEobj@F_FMSY[, 
            mm, ]), na.rm = T) * 100, 1)
        P10[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
        P50[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 0.5, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
        P100[mm] <- round(sum(MSEobj@B_BMSY[, mm, ] <= 1, na.rm = T)/prod(dim(MSEobj@B_BMSY[, 
            mm, ])) * 100, 1)
    }
    coly <- rep(brewer.pal(8,"Dark2"), 20)##[1:length(MSEobj@MPs[1:MSEobj@nMPs])]    
    coly[MSEobj@MPs[1:MSEobj@nMPs] %in% c("AvC", "curE", "FMSYref")] <- "black"    
##    old_par <- par(mfrow = c(2, 2), mai = c(0.9, 1, 0.1, 0.1), 
##        omi = c(0.1, 0.1, 0.4, 0))
##    tradeoffplotTKM(POF, Yd, "Prob. of overfishing (%)", "Relative yield", 
##                 MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
##    tradeoffplotTKM(P100, Yd, "Prob. biomass < BMSY (%)", "Relative yield", 
##        MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
    tradeoffplotTKM(P50, Yd, "Prob. biomass < 0.5BMSY (%)", "Relative yield", 
                    MSEobj@MPs[1:MSEobj@nMPs], vl = c(5,50), hl = 100)
##    tradeoffplotTKM(P10, Yd, "Prob. biomass < 0.1BMSY (%)", "Relative yield", 
##                 MSEobj@MPs[1:MSEobj@nMPs], vl = 50, hl = 100)
##    legend("bottomright", legend=MSEobj@MPs[1:MSEobj@nMPs], col = coly, bty="n",pch=16, cex=0.8)    

        ##        mtext(deparse(substitute(MSEobj)), 3, outer = T, line = 0.3, font = 2)
        mtext(main, 3, outer = F, line = 0.3, font = 2)
##    if (!is.na(nam) & !is.character(nam)) 
##        mtext(MSEobj@Name, 3, outer = T, line = 0.3, font = 2)
##    if (!is.na(nam) & is.character(nam)) 
##        mtext(nam, 3, outer = T, line = 0.3, font = 2)
##    par(old_par)
}
