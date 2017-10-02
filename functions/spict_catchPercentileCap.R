#' @name spict_catchPercentileCap
#' @title SPiCT assessment with specified percentiles from the predicted catch distribution
#' 
#' @details SPiCT assessment is done using catch and relative biomass index observations. 
#' Stock status estimates are used to set the TAC for the next year, equal to the
#' catch that corresponds to fishing mortality equal to 80\% of Fmsy.
#'
#' @param x A position in a data-limited mehods data object
#' @param DLM_data A data-limited methods data object (see DLMtool)
#' @param reps The number of stochastic samples of the TAC recommendation
#'
#' @return A numeric vector of TAC recommendations
#' @export
#'
#' @examples
#' \dontrun{
#' library(DLMtool)
#' 
#' ## Put together an operating model from the available DLM toolkit examples
#' stock <- Herring
#' Fleet.example <- Generic_IncE
#' Observation.example <- Precise_Unbiased
#' 
#' ## Remove changes in life history parameters
#' stock@Mgrad <- c(0,0)
#' stock@Kgrad <- c(0,0)
#' stock@Linfgrad <- c(0,0)
#' stock@Prob_staying <- c(1,1)
#' 
#' ## Set the depletion level 
#' stock@D <- c(0.3, 0.4)
#'
#' OM.example <- new("OM", Stock = stock, Fleet = Fleet.example, 
#'                   Obs = Observation.example)
#'
#' MP.vec <- c("spict_catchPercentile50")
#' 
#' MSE.example <- runMSE(OM.example, MPs = MP.vec,
#'                       interval = 1, reps = 100, timelimit = 150, CheckMPs = FALSE)
#' }


spict_catchPercentile05cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.05, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)


spict_catchPercentile10cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.1, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile15cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.15, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile20cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.2, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile25cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.25, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile30cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.3, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile35cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.35, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile40cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.4, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile45cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.45, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)



spict_catchPercentile50cap <- structure(
    function(x, Data, reps = 1, percentileC = 0.5, percentileFmsy=NA, cap=TRUE){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 16,
                    do.sd.report=TRUE,
                    getReportCovariance = FALSE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error") || rep$opt$convergence != 0) {
            TAC <- rep(NA, reps)
        } else {
            if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-percentileFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm
            } else {
                red <- 1
            }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                if(cap){
                    idx <- rep$inp$indpred[1]                    
                    blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                    bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)
