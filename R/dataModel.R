#' Create a DataModel object, describing the relationship of your data to the SEIR model
#' 
#' @param Y  A matrix with T rows and n columns, for time points and spatial locations respectively.
#' @param type  A string equal to "identity", "overdispersion", or "fractional"
#' @param compartment  A string equal to "I_star" or "R_star" 
#' @param cumulative  A logical value indicating whether the data is reported on a cumulative scale
#' @param params  A list of optional parameters 
#'         
#' @return an object of type \code{\link{DataModel}} 
#'
#' @details 
#'  A fundamental task in building hierarchical models is to describe the way in which
#'  the observed data relates to the underlying model. Here, we assume that the 
#'  matrix of observed values, \code{Y}, is related to one of two epidemic compartments:
#'  I_star or R_star. In the notation of Brown, Porter, and Oleson 2015, these correspond
#'  to the $I^*$ and $R^*$ compartments, which catalog the newly infectious and recovered
#'  individuals, respectively. For a discussion of these definitions, in addition to a more
#'  detailed description of the overall spatial SEIR framework, please refer to that work. 
#'  Depending on the data model type, additional parameters may be required:
#'  \itemize{
#'  \item{phi: }{An overdispersion parameter, required for the 'overdispersion' data model type.}
#'  \item{report_fraction}{The (scalar) estimated reporting fraction for the epidemic. 
#'  Required for the 'fractional' data model.}
#'  \item{report_fraction_ess}{The effective sample size associated with the required reporting fraction
#'   in the 'fractional' data model.}}
#'
#'
#' @examples dataModel = DataModel(rpois(100, 10)) 
#' @export
DataModel = function(Y, type = c("identity", "overdispersion"), 
                     compartment = c("I_star", "R_star", "I"), 
                     cumulative=FALSE,
                     params=list())
{
    type = type[1] 
    compartment = compartment[1] 
    checkArgument("params", 
                  mustHaveClass("list"),
                  validateIf(type == "overdispersion", 
                             mustHaveMember("phi")),
                  validateIf(type == "fractional",
                             mustHaveMember("report_fraction")),
                  validateIf(type == "fractional",
                             mustHaveMember("report_fraction_ess")) 
                  )
    if (class(Y) != "matrix")
    {
        Y = as.matrix(Y)
    }


    phi <- ifelse("phi" %in% names(params), params$phi, -1)
    report_fraction <- ifelse("report_fraction" %in% names(params), params$report_fraction, -1)
    report_fraction_ess <- ifelse("report_fraction_ess" %in% names(params), params$report_fraction_ess, -1)
    weights <- ifelse("weights" %in% names(params), params$weights, rep(1, ncol(Y)))

    checkArgument("report_fraction",
                  validateIf(type=="fractional",
                             mustBeInRange(0,1,FALSE, TRUE)
                  ))
    checkArgument("report_fraction_ess",
                  validateIf(type=="fractional",
                             mustBeInRange(lower=1)
                 ))


    if (type == "overdispersion" && (!("phi" %in% names(params)) 
                                     || class(params[["phi"]]) != "numeric"))
    {
        stop("Numeric phi required for overdispersion models")
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
    else if (length(cumulative) != 1 || class(cumulative) != "logical")
    {
        stop("The cumulative argument must be a logical value of length 1.")
    }
    
    if (length(weights) != ncol(Y) || any(!is.finite(weights))){
        stop("weights, if specified, must be of length ncol(Y) and finite.")
    }
    na_mask = is.na(Y)
    Y[na_mask] = -Inf
    structure(list("Y"=Y, 
                   "type"=type,
                   "compartment"=compartment,
                   "cumulative"=cumulative,
                   "phi"=phi,
                   "report_fraction"=report_fraction,
                   "report_fraction_ess"=report_fraction_ess,
                   "na_mask" = na_mask,
                   "weights" = weights), class = "DataModel")
}

