#' Create a DistanceModel object describing the network structure of the 
#' population under study. The TDistanceModel function accepts lagged
#` contact.

#' @param distanceList a list of square, symmetric distance matrices.
#' @param laggedDistanceList a list of lists, of dimension T by k, where
#'  T is the number of time points and k is the number of time points being
#'  carried forward.
#' @param scaleMode an optional argument specifying the type of 
#' preprocessing needed.
#' @param priorAlpha the first shape parameter for the beta distributed
#' autocorrelation terms
#' @param priorBeta the second shape parameter for the beta distributed
#' autocorrelation terms
#' @return an object of type \code{\link{TDistanceModel}}
#' @details
#'  In stochastic spatial SEIR models as specified in Brown et al. 2015, 
#'  populations are divided into homogeneous groups, or locations, with
#'  heterogeneous mixing between groups. This is accomplished 
#'  using a distance matrix parameterization, in which some number of 
#'  square, symmetric distance matrices are constructed, each of which 
#'  receives a spatial autocorrelation parameter. 
#' 
#'  Care must be taken to specify reasonable prior parameters for such
#'  terms, as well as in the construction and scaling of the distance 
#'  matrices; it is certainly possible to construct an overspecified
#'  model, and to correspondingly bias inference about other important
#'  exposure process terms. 
#' 
#' @examples distanceModel <- DistanceModel(list(1-diag(4)))
#' @export
TDistanceModel = function(distanceList, 
                         laggedDistanceList,
                         scaleMode = c("none","rowscale","invsqrt"),
                         priorAlpha=1.0,
                         priorBeta=1.0)
{
    scaleMode = scaleMode[1]
    rowScale = function(mat)
    {
        mat/matrix(apply(mat,1,sum), nrow = nrow(mat), ncol = ncol(mat))
    }

    invSqrt = function(mat)
    {
        matrix(ifelse(mat == 0, 0, 1/sqrt(mat)), nrow = nrow(mat), ncol = ncol(mat))
    }


    # Check if we've got a single matrix
    if (class(distanceList) == "matrix")
    {
        distanceList = list(distanceList)
    }
    # Make sure that at this point we have a list of matrices. 

    if (class(distanceList) != "list")
    {
        stop("Error: distanceList must be a list of matrices.")
    }

    if (class(laggedDistanceList) != "list")
    {
        stop("Error: laggedDistanceList must be a list of matrices.")
    }
    nLags <- ifelse(length(laggedDistanceList) == 0, 0, 
                    length(laggedDistanceList[[1]]))


    # Check for valid matrices
    distanceDim = NA
    for (i in 1:length(distanceList))
    {
        if (class(distanceList[[i]]) != "matrix")
        {
            stop("Distance metrics must be matrices.")
        }
        distanceDim = dim(distanceList[[1]])
        newDim = dim(distanceList[[i]])
        if (any(newDim != distanceDim) || (newDim[1] != newDim[2]))
        {
            stop("Distance matrices must be square and of the same dimension.")
        }
        if (any(is.na(distanceList[[i]])))
        {
            stop("NA's are not allowed in distance matrices.")
        }
        if (scaleMode == "rowscale")
        {
            distanceList[[i]] = rowScale(distanceList[[i]])
        }
        if (scaleMode == "invsqrt")
        {
            distanceList[[i]] = invSqrt(distanceList[[i]])
        }
    } 

    structure(list("distanceList" = distanceList,
                   "laggedDistanceList" = laggedDistanceList,
                   "len" = nLags*length(laggedDistanceList) + length(distanceList),
                   "priorAlpha" = priorAlpha,
                   "priorBeta" = priorBeta), class = "DistanceModel")
}


