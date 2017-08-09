#' Print a \code{summary.SpatialSEIRModel} object
#'
#' @param x a \code{summary.SpatialSEIRModel} object
#' @param \dots not used
#' @export
print.summary.SpatialSEIRModel = function(x, ...)
{
    nl = "\n\n"
    if (x$transitionParams > 1){
        cat(paste("Summary: SEIR Model", nl)) 
    }
    else{
        cat(paste("Summary: SIR Model", nl)) 
    }
    cat(paste("Locations: ", x$nLoc, "\n",
              "Time Points: ", x$nTpt,"\n", sep = ""))
    cat(paste("Data Model Parameters: ", 
              x$dataModelParams, "\n", sep = ""))
    cat(paste("Exposure Process Parameters: ", 
              x$exposureParams, "\n", sep = ""))
    cat(paste("Reinfection Model Parameters: ", 
             (x$reinfectionParams), "\n", sep = ""))
    cat(paste("Spatial Parameters: ", 
              (x$spatialParams), "\n", sep = ""))
    cat(paste("Transition Parameters: ", 
              (x$transitionParams), "\n", sep = ""))

    cat(nl)
    cat("Parameter Estimates:\n")
    print(round(x$parameterEstimates, 3))
    cat("\n") 
}
