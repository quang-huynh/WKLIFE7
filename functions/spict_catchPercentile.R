#' @name spict_catchPercentile
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

spict_catchPercentile05 <- structure(
    function(x, Data, reps = 1, quant = 0.05){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)


spict_catchPercentile10 <- structure(
    function(x, Data, reps = 1, quant = 0.1){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile15 <- structure(
    function(x, Data, reps = 1, quant = 0.15){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile20 <- structure(
    function(x, Data, reps = 1, quant = 0.2){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile25 <- structure(
    function(x, Data, reps = 1, quant = 0.25){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()        
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile30 <- structure(
    function(x, Data, reps = 1, quant = 0.3){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc();        
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile35 <- structure(
    function(x, Data, reps = 1, quant = 0.35){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile40 <- structure(
    function(x, Data, reps = 1, quant = 0.40){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile45 <- structure(
    function(x, Data, reps = 1, quant = 0.45){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile50 <- structure(
    function(x, Data, reps = 1, quant = 0.50){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile55 <- structure(
    function(x, Data, reps = 1, quant = 0.55){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile60 <- structure(
    function(x, Data, reps = 1, quant = 0.6){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile65 <- structure(
    function(x, Data, reps = 1, quant = 0.65){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile70 <- structure(
    function(x, Data, reps = 1, quant = 0.70){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile75 <- structure(
    function(x, Data, reps = 1, quant = 0.75){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile80 <- structure(
    function(x, Data, reps = 1, quant = 0.8){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile85 <- structure(
    function(x, Data, reps = 1, quant = 0.85){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile90 <- structure(
    function(x, Data, reps = 1, quant = 0.9){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)

spict_catchPercentile95 <- structure(
    function(x, Data, reps = 1, quant = 0.95){
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
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }

        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        c(tacTemp, rep(NA, reps-1))
        
    },
    class = "Output"
)
