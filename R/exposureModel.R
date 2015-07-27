# exposureModel module helper function
ExposureModel = function(X,nTpt, nLoc,betaPriorPrecision=NA,
                              betaPriorMean=NA,offset=NA)
{
    nBeta = ncol(X)
    if (class(X) != "matrix")
    {
        print("Warning: X should be a matrix.")
    }
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
    structure(list(X=X,
                   nTpt=nTpt,
                   nLoc=nLoc,
                   betaPriorMean=betaPriorMean,
                   betaPriorPrecision=betaPriorPrecision,
                   offset=offset), class="ExposureModel")
}


