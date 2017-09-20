#' @name SPiCT_catchQuant
#' @title SPiCT assessment with specified quantiles from the predicted catch distribution
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
#' MP.vec <- c("SPiCT_catchQuant")
#' 
#' MSE.example <- runMSE(OM.example, MPs = MP.vec, nsim = 200, proyears = 20,
#'                       interval = 1, reps = 100, timelimit = 150, CheckMPs = FALSE)
#' }

SPiCT_catchQuant <- structure(
    function(x, Data, reps = 1, quant = 0.5){
        dependencies <- "Data@Year, Data@Cat, Data@Ind"
        time <- Data@Year
        Catch <- Data@Cat[x,]
        Index <- Data@Ind[x,]
        inp <- list(timeC=time, obsC=Catch, 
                    timeI=time, obsI=Index,
                    ## timepredc = max(time) + 1,
                    dteuler = 1 / 4,
                    do.sd.report=TRUE,
                    getReportCovariance = TRUE)
        rep <- try(spict::fit.spict(inp))
        if(is(rep, "try-error")) {
            TAC <- rep(NA, reps)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, reps)
            } else {
                TAC <- exp(qnorm(quant, predcatch[2], predcatch[4]))
            }
        }
        DLMtool:::TACfilter(TAC)
    },
    class = "Output"
)
