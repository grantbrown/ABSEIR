#' Produce a summary of the results of fitting a SEIR model
#' @param object a \code{\link{SpatialSEIRModel}} object
#' @param \dots not used 
#' @return a summary.SpatialSEIRModel object
#' @export
summary.SpatialSEIRModel = function(object, ...)
{
    nLoc = ncol(object$modelComponents$data_model$Y)
    nTpt = nrow(object$modelComponents$data_model$Y)

    nLags = length(object$modelComponents$distance_model$laggedDistanceList[[1]]) 
    transitionMode = object$modelComponents$transition_priors$mode
    hasSpatial = (object$modelComponents$exposure_model$nLoc > 1) 
    hasReinfection = 
        (object$modelComponents$reinfection_model$integerMode != 3) 
    hasReportFraction = object$modelComponents$data_model$type == "fractional"
    hasLatent = (object$modelComponents$transition_priors$enable_latent)

    qtiles = t(apply(object$param.samples, 2, quantile, 
                   probs = c(0.025, 0.975)))

    means = apply(object$param.samples, 2, mean)
    sds = apply(object$param.samples, 2, sd)

    outMatrix = cbind(means, sds, qtiles)
    rownames(outMatrix) = colnames(object$param.samples)
    colnames(outMatrix) = c("Mean", "SD", "95% LB", "95% UB")
    structure(list(modelComponents=object$modelComponents,
       parameterEstimates=outMatrix,
       nLoc = nLoc,
       nTpt = nTpt,
       exposureParams = ncol(object$modelComponents$exposure_model$X),
       reinfectionParams = Ifelse(hasReinfection, 
            ncol(object$modelComponents$reinfection_model$X_prs), 
            0),
       spatialParams = Ifelse(hasSpatial, 
            length(object$modelComponents$distance_model$distanceList) + nLags,
            0),
       transitionParams = Ifelse(transitionMode == "exponential", 
                                 (ifelse(hasLatent, 2, 1)),
                          Ifelse(transitionMode == "weibull", 4, 0)),
       dataModelParams = Ifelse(hasReportFraction, 1, 0)
       ), class = "summary.SpatialSEIRModel")
}

