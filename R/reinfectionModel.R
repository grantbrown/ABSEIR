# reinfectionModel module helper function
ReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, 
                                 betaPrs=NA , priorPrecision = NA, priorMean = NA)
{
    integerMode = ifelse(reinfectMode[1] == "SEIR", 3, 
                  ifelse(reinfectMode[1] == "SEIRS", 1, 
                  ifelse(reinfectMode[1] == "Fixed", 2, NA)))
    if (is.na(integerMode))
    {
        stop(paste("Invalid mode: ", reinfectMode[1], sep = ""))
    }
    if (integerMode != 3 && (is.na(X_prs) || is.na(betaPrs)))
    {
        stop("If reinfection mode is not SEIR, X_prs and betaPrs must be supplied.")
    }
    if (integerMode == 1 && (is.na(priorPrecision) || is.na(priorMean)))
    {
        stop("If reinfection parameters are going to be estimated, priorPrecision and priorMean must be specified.")
    }

    if (integerMode != 3)
    {
        if (class(X_prs) != "matrix")
        {
            X_prs = as.matrix(X_prs)
        }
        if (all(is.na(priorPrecision)))
        {
           priorPrecision = rep(0.1, ncol(X_prs))     
        }
        else if (length(priorPrecision) == 1)
        {
            priorPrecision = rep(priorPrecision, ncol(X_prs))               
        }
        if (all(is.na(priorMean)))
        {
            priorMean = rep(0, ncol(X_prs))
        }
        else if (length(priorMean) == 1)
        {
            priorMean = rep(priorMean, ncol(X_prs))
        }
    }

    structure(list("integerMode" = integerMode,
                   "X_prs" = X_prs,
                   "priorPrecision" = priorPrecision,
                   "priorMean" = priorMean,
                   "betaPrs" = betaPrs
                   ), class = "ReinfectionModel")
}


