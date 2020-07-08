#' Create an ExposureModel object, describing an epidemic intensity process
#'
#' @param X an $(n*T)$ by $p$ design matrix, where $n$ is the number
#' of locations, $T$ is the number of time points, and $p$ is the number
#' of exposure process parameters. Each column corresponds to a parameter,
#' while each block of $T$ rows corresponds to the time series of covariate
#' values associated with a location. 
#' @param nTpt the number of time points, $T$
#' @param nLoc the number of locations, $N$
#' @param betaPriorPrecision the prior precisions of the $p$ exposure 
#' process parameters
#' @param betaPriorMean the prior means of the $p$ exposure process
#' parameters
#' @param offset a vector of $T$ temporal offset terms, capturing the 
#' relative ammount of aggregated continuous time corresponding to each
#' recorded discrete time point. 
#'
#' @details
#' The exposure process allows the inclusion of both spatially and temporally
#' varying covariates, as well as invariant quantities (intercepts, demographic 
#' features etc.).  
#' @examples exposureModel <- ExposureModel(cbind(1, (1:25)/25),
#'                                          nTpt = 25, nLoc = 1,
#'                                          betaPriorPrecision = 0.1,
#'                                          betaPriorMean = 0)
#' @export
ExposureModel = function(X,nTpt, nLoc,betaPriorPrecision=NA,
                              betaPriorMean=NA,offset=NA)
{
    nBeta = ncol(X)
    checkArgument("X", mustHaveClass("matrix"))

    if (nrow(X) != nTpt * nLoc){
        stop("Invalid data dimensions: X should be a matrix composed of nLoc row-wise blocks of dimension nTpt*p.") 
    } 
    if (length(betaPriorPrecision) == 1 && is.na(betaPriorPrecision))
    {
        print("No prior precision specified, using 0.1")
        betaPriorPrecision = rep(0.1, nBeta)
    }
    else if (length(betaPriorPrecision) == 1)
    {
        betaPriorPrecision = rep(betaPriorPrecision, nBeta)
    }
    if (length(betaPriorMean) == 1 && is.na(betaPriorMean))
    {
        print("No prior mean specified, using zero.")
        betaPriorMean = rep(0, nBeta)
    }
    else if (length(betaPriorMean) == 1)
    {
        betaPriorMean = rep(betaPriorMean, nBeta)
    }
    if (length(offset) == 1 && is.na(offset))
    {
        print("Assuming equally spaced count data.") 
    }
    else if (length(offset) == 1){
        offset = rep(offset, nTpt)
    }
    else if (length(offset) != nTpt){
        stop("Offset must be the of length nTpt")
    }
    if (length(offset) != 1 && max(abs(offset - as.integer(offset))) > 0)
    {
        stop("Offsets must be integers - we need a mapping onto discrete time units.")
    }
    structure(list(X=X,
                   nTpt=nTpt,
                   nLoc=nLoc,
                   betaPriorMean=betaPriorMean,
                   betaPriorPrecision=betaPriorPrecision,
                   offset=offset), class="ExposureModel")
}


