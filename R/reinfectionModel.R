#' Create a ReinfectionModel object, which determines whether and how individuals
#' may return to the susceptible population after an infection runs its course. 
#' 
#' @param reinfectMode a string equal to one of: \code{"SEIR"}
#' , \code{"SEIRS"}, or \code{"Fixed"}, SEIR models include no
#' reinfection. SEIRS models include a reinfection process,
#' and Fixed models allow reinfection, but do not estimate the parameters
#' involved. Currently, only SEIR and SEIRS modes are supported. 
#' @param X_prs a $T$ by $p$ design matrix, which must be specified for SEIRS and
#' Fixed models. 
#' @param priorPrecision the prior precisions for each of the reinfection parameters
#' @param priorMean the prior means for each of the reinfection parameters
#' @details
#'  SEIR models are probably the most commonly employed version of compartmental
#'  epidemic models today, and are specified using a simple SEIR reinfection model
#'  (i.e., a model which introduces no reinfection process). On the other hand, 
#'  there are some diseases which must be modeled more flexibly, for individuals
#'  are not permanently (over the course of a typical study period) removed from
#'  the susceptible population by immunity or death after an infection. 
#' @examples reinfectionModel = ReinfectionModel("SEIR")
#' 

ReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, 
                                 priorPrecision = NA, priorMean = NA)
{
    integerMode = ifelse(reinfectMode[1] == "SEIR", 3, 
                  ifelse(reinfectMode[1] == "SEIRS", 1, 
                  ifelse(reinfectMode[1] == "Fixed", 2, NA)))
    if (is.na(integerMode))
    {
        stop(paste("Invalid mode: ", reinfectMode[1], sep = ""))
    }
    if (integerMode != 3 && (is.na(X_prs)))
    {
        stop("If reinfection mode is not SEIR, X_prs must be supplied.")
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
                   "priorMean" = priorMean
                   ), class = "ReinfectionModel")
}


