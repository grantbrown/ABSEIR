# SEIRModel module helper function
SpatialSEIRModel = function(data_model,
                          exposure_model,
                          reinfection_model,
                          distance_model,
                          transition_priors,
                          initial_value_container,
                          sampling_control,
                          samples=100,
                          verbose=FALSE)
{

    if (class(data_model) != "DataModel")
    {
        stop(paste("Expected: DataModel. Received: ", 
                   class(data_model)))
    }
    if (class(exposure_model) != "ExposureModel")
    {
        stop(paste("Expected: ExposureModel Received: ", 
                   class(exposure_model)))
    }
    if (class(reinfection_model) != "ReinfectionModel")
    {
        stop(paste("Expected: ReinfectionModel Received: ", 
                   class(reinfection_model)))
    }
    if (class(distance_model) != "DistanceModel")
    {
        stop(paste("Expected: DistanceModel Received: ", 
                   class(distance_model)))
    }
    if (class(transition_priors) != "TransitionPriors")
    {
        stop(paste("Expected: TransitionPriors Received: ", 
                   class(transition_priors)))
    }
    if (class(sampling_control) != "SamplingControl")
    {
        stop(paste("Expected: SamplingControl Received: ", 
                   class(sampling_control)))
    }

    modelResults = list()
    modelComponents = list()
    if (verbose){cat("Initializing Model Components\n")}
    hasSpatial = (ncol(data_model$Y) > 1) 
    hasReinfection = (reinfection_model$integerMode != 3) 
    result = tryCatch({
        if (verbose) cat("...Building data model\n")
        modelComponents[["dataModel"]] = new(dataModel, data_model$Y,
                                            data_model$type,
                                            data_model$compartment,
                                            data_model$phi,
                                            data_model$na_mask*1)

        if (verbose) cat("...Building distance model\n")
        modelComponents[["distanceModel"]] = new(distanceModel)
        for (i in 1:length(distance_model$distanceList))
        {
            modelComponents[["distanceModel"]]$addDistanceMatrix(
                distance_model$distanceList[[i]]
            )
        }
        modelComponents[["distanceModel"]]$setPriorParameters(
            distance_model$priorAlpha,
            distance_model$priorBeta
        )

        if (verbose) cat("...Building exposure model\n")
        modelComponents[["exposureModel"]] = new(
            exposureModel, 
            exposure_model$X,
            exposure_model$nTpt,
            exposure_model$nLoc,
            exposure_model$betaPriorMean,
            exposure_model$betaPriorPrecision
        )
        if (!all(is.na(exposure_model$offset)))
        {
            modelComponents[["exposureModel"]]$offsets = (
                exposure_model$offset
            )
        }

        if (verbose) cat("...Building initial value container\n")
        modelComponents[["initialValueContainer"]] = new(initialValueContainer)
        modelComponents[["initialValueContainer"]]$setInitialValues(
            initial_value_container$S0,
            initial_value_container$E0,
            initial_value_container$I0,
            initial_value_container$R0
        )

        if (verbose) cat("...Building reinfection model\n") 
        modelComponents[["reinfectionModel"]] = new(
            reinfectionModel, 
            reinfection_model$integerMode
        )
        if (reinfection_model$integerMode != 3)
        {
            modelComponents[["reinfectionModel"]]$buildReinfectionModel(
                reinfection_model$X_prs, 
                reinfection_model$priorMean, 
                reinfection_model$priorPrecision
        );
        }

        if (verbose) cat("...Building sampling control model\n") 
        modelComponents[["samplingControl"]] = new (
            samplingControl, 
            c(sampling_control$sim_width, sampling_control$seed,
              sampling_control$n_cores,sampling_control$algorithm, 
              sampling_control$batch_size,sampling_control$epochs, 
              sampling_control$max_batches),
            c(sampling_control$acceptance_fraction, sampling_control$shrinkage)
        )

        if (verbose) cat("...Building transition priors\n") 
        modelComponents[["transitionPriors"]] = new(transitionPriors)
        modelComponents[["transitionPriors"]]$setPriorsFromProbabilities(
            transition_priors$p_ei,
            transition_priors$p_ir,
            transition_priors$p_ei_ess,
            transition_priors$p_ir_ess
        )

        if (verbose) cat("...Preparing model object\n")
        modelComponents[["SEIR_model"]] = new( 
            spatialSEIRModel, 
            modelComponents[["dataModel"]],
            modelComponents[["exposureModel"]],
            modelComponents[["reinfectionModel"]],
            modelComponents[["distanceModel"]],
            modelComponents[["transitionPriors"]],
            modelComponents[["initialValueContainer"]],
            modelComponents[["samplingControl"]]
        )
        if (verbose) cat("Running main simulation\n")
        rslt = modelComponents[["SEIR_model"]]$sample(samples, verbose*1)
        if (verbose) cat("Simulation complete\n")

        epsilon = rslt$result
        params = rslt$params

        cnames = paste("Beta_SE_", 
                       1:ncol(exposure_model$X), sep = "")
        if (hasReinfection)
        {
            cnames = c(cnames, paste("Beta_RS_", 
                                     1:ncol(reinfection_model$X_prs), 
                                     sep = "")
            )
        }
        if (hasSpatial)
        {
            cnames = c(cnames, 
                       paste("rho_", 
                            1:length(distance_model$distanceList), 
                             sep = "")
            )
        }
        cnames = c(cnames, "gamma_EI", "gamma_IR")
        colnames(params) = cnames 
        
        modelResults[["param.samples"]] = params
        modelResults[["epsilon"]] = epsilon
        modelResults[["modelComponents"]] = list(
                     data_model = data_model,
                     exposure_model = exposure_model,
                     reinfection_model = reinfection_model,
                     distance_model = distance_model,
                     transition_priors = transition_priors,
                     initial_value_container = initial_value_container,
                     sampling_control = sampling_control
        ) 
    },
    warning=function(w)
    {
        cat("Warnings generated:\n")
        print(w)
    },
    error=function(e)
    {
        cat("Errors generated:\n")
        print(e)
        stop("Aborting")
    },
    finally = { 
        if (exists("modelComponents"))
        {
            rm(modelComponents)
        }
    })
    return(structure(modelResults, class = "SpatialSEIRModel"))
}

#' Plot a graphical summary of the marginal posterior distribution of SEIR model parameters
#' @param x a \code{\link{SpatialSEIRModel}} object
#' @param \dots additional arguments to be passed to plotting functions. 
#' @examples \dontrun{plot(modelObject)}
plot.SpatialSEIRModel = function(x, ...)
{
    for (i in 1:ncol(x$param.samples))
    {
        cat ("Press return to see next plot:")
        hist(x$param.samples[,i], 
             main = colnames(x$param.samples)[i],
             xlab = "Parameter Value",
             ylab = "MPD",
             freq=FALSE,
             ...
             )
        if (i != ncol(x$param.samples)) line <- readline()
    }


}

#' Produce a summary of the results of fitting a SEIR model
#' @param object a \code{\link{SpatialSEIRModel}} object
#' @param \dots not used 
#' @return a summary.SpatialSEIRModel object
summary.SpatialSEIRModel = function(object, ...)
{
    nLoc = ncol(object$modelComponents$data_model$Y)
    nTpt = nrow(object$modelComponents$data_model$Y)

    hasSpatial = (object$modelComponents$exposure_model$nLoc > 1) 
    hasReinfection = 
        (object$modelComponents$reinfection_model$integerMode != 3) 

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
       reinfectionParams = ifelse(hasReinfection, 
            ncol(object$modelComponents$reinfection_model$X_prs), 
            0),
       spatialParams = ifelse(hasSpatial, 
            length(object$modelComponents$distance_model$distanceList),
            0)
       ), class = "summary.SpatialSEIRModel")
}

#' Print a \code{summary.SpatialSEIRModel} object
#' @param x a \code{summary.SpatialSEIRModel} object
#' @param \dots not used
print.summary.SpatialSEIRModel = function(x, ...)
{
    nl = "\n\n"
    cat(paste("Summary: SEIR Model", nl)) 
    cat(paste("Locations: ", x$nLoc, "\n",
              "Time Points: ", x$nTpt,"\n", sep = ""))
    cat(paste("Exposure Process Parameters: ", 
              x$exposureParams, "\n", sep = ""))
    cat(paste("Reinfection Model Parameters: ", 
             (x$reinfectionParams), "\n", sep = ""))
    cat(paste("Spatial Parameters: ", 
              (x$spatialParams), "\n", sep = ""))
    cat(nl)
    cat("Parameter Estimates:\n")
    print(round(x$parameterEstimates, 3))
    cat("\n") 
}
