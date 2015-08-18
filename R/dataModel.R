#' Create a DataModel object, describing the relationship of your data to the SEIR model
#' 
#' @param Y a matrix with T rows and n columns, for time points and spatial locations respectively.
#' @param type a string equal to "identity" or "overdispersion"
#' @param compartment a string equal to "I_star" or "R_star" 
#' @param phi optional overdispersion parameter
#' @return an object of type \code{\link{DataModel}} 
#' @details
#'  A fundamental task in building hierarchical models is to describe the way in which
#'  the observed data relates to the underlying model. Here, we assume that the 
#'  matrix of observed values, \code{Y}, is related to one of two epidemic compartments:
#'  I_star or R_star. In the notation of Brown, Porter, and Oleson 2015, these correspond
#'  to the $I^*$ and $R^*$ compartments, which catalog the newly infectious and recovered
#'  individuals, respectively. For a discussion of these definitions, in addition to a more
#'  detailed description of the overall spatial SEIR framework, please refer to that work. 
#' @examples print("Examples go here") 

DataModel = function(Y, type = c("identity", "overdispersion"), compartment = c("I_star", "R_star", "I"), phi = NA)
{
    type = type[1] 
    compartment = compartment[1] 
    if (compartment != "I_star")
    {
        stop(paste("Currently, only the I_star compartment is supported.", 
             " Alternate models are planned, so if you find yourself",
             " in need of such features, please email grant-brown@uiowa.edu"))
    }
    if (class(Y) != "matrix")
    {
        Y = as.matrix(Y)
    }
    if (length(phi) == 1 && is.na(phi) && type != "identity")
    {
        stop("Must specify overdispersion parameter (phi) for non-identity data model.")
    }
    else if (length(phi) == 1 && is.na(phi))
    {
        phi = -1
    }
    else if (class(phi) != "numeric" || length(phi) != 1)
    {
        stop("Non identy data model currently requires a single overdispersion parameter") 
    }
    na_mask = is.na(Y)
    Y[na_mask] = -Inf
    structure(list("Y"=Y, 
                   "type"=type,
                   "compartment"=compartment,
                   "phi"=phi,
                   "na_mask" = na_mask), class = "DataModel")
}

