# SEIRModel module helper function
SpatialSEIRModel = function(dataModelInstance,
                          exposureModelInstance,
                          reinfectionModelInstance,
                          distanceModelInstance,
                          transitionPriorsInstance,
                          initialValueContainerInstance,
                          samplingControlInstance,
                          samples=100,
                          acceptFraction=1e-4,
                          batchSize=50000,
                          verbose=FALSE)
{

    if (class(dataModelInstance) != "DataModel")
    {
        stop(paste("Expected: DataModel. Received: ", 
                   class(dataModelInstance)))
    }
    if (class(exposureModelInstance) != "ExposureModel")
    {
        stop(paste("Expected: ExposureModel Received: ", 
                   class(exposureModelInstance)))
    }
    if (class(reinfectionModelInstance) != "ReinfectionModel")
    {
        stop(paste("Expected: ReinfectionModel Received: ", 
                   class(reinfectionModelInstance)))
    }
    if (class(distanceModelInstance) != "DistanceModel")
    {
        stop(paste("Expected: DistanceModel Received: ", 
                   class(distanceModelInstance)))
    }
    if (class(transitionPriorsInstance) != "TransitionPriors")
    {
        stop(paste("Expected: TransitionPriors Received: ", 
                   class(transitionPriorsInstance)))
    }
    if (class(samplingControlInstance) != "SamplingControl")
    {
        stop(paste("Expected: SamplingControl Received: ", 
                   class(samplingControlInstance)))
    }

    modelResults = list()
    modelComponents = list()
    if (verbose){cat("Initializing Model Components\n")}
    hasSpatial = (ncol(dataModelInstance$Y) > 1) 
    hasReinfection = (reinfectionModelInstance$integerMode != 3) 
    result = tryCatch({
        if (verbose) cat("...Building data model\n")
        modelComponents[["dataModel"]] = new(dataModel, dataModelInstance$Y,
                                            dataModelInstance$type,
                                            dataModelInstance$compartment,
                                            dataModelInstance$phi)

        if (verbose) cat("...Building distance model\n")
        modelComponents[["distanceModel"]] = new(distanceModel)
        for (i in 1:length(distanceModelInstance$distanceList))
        {
            modelComponents[["distanceModel"]]$addDistanceMatrix(
                distanceModelInstance$distanceList[[i]]
            )
        }
        modelComponents[["distanceModel"]]$setPriorParameters(
            distanceModelInstance$priorAlpha,
            distanceModelInstance$priorBeta
        )

        if (verbose) cat("...Building exposure model\n")
        modelComponents[["exposureModel"]] = new(
            exposureModel, 
            exposureModelInstance$X,
            exposureModelInstance$nTpt,
            exposureModelInstance$nLoc,
            exposureModelInstance$betaPriorMean,
            exposureModelInstance$betaPriorPrecision
        )
        if (!all(is.na(exposureModelInstance$offset)))
        {
            modelComponents[["exposureModel"]]$offsets = (
                exposureModelInstance$offset
            )
        }

        if (verbose) cat("...Building initial value container\n")
        modelComponents[["initialValueContainer"]] = new(initialValueContainer)
        modelComponents[["initialValueContainer"]]$setInitialValues(
            initialValueContainerInstance$S0,
            initialValueContainerInstance$E0,
            initialValueContainerInstance$I0,
            initialValueContainerInstance$R0
        )

        if (verbose) cat("...Building reinfection model\n") 
        modelComponents[["reinfectionModel"]] = new(
            reinfectionModel, 
            reinfectionModelInstance$integerMode
        )
        if (reinfectionModelInstance$integerMode != 3)
        {
            modelComponents[["reinfectionModel"]]$buildReinfectionModel(
                reinfectionModelInstance$X_prs, 
                reinfectionModelInstance$priorMean, 
                reinfectionModelInstance$priorPrecision
        );
        }

        if (verbose) cat("...Building sampling control model\n") 
        modelComponents[["samplingControl"]] = new (
            samplingControl, 
            samplingControlInstance$sim_width,
            samplingControlInstance$seed,
            samplingControlInstance$n_cores
        )

        if (verbose) cat("...Building transition priors\n") 
        modelComponents[["transitionPriors"]] = new(transitionPriors)
        modelComponents[["transitionPriors"]]$setPriorsFromProbabilities(
            transitionPriorsInstance$p_ei,
            transitionPriorsInstance$p_ir,
            transitionPriorsInstance$p_ei_ess,
            transitionPriorsInstance$p_ir_ess
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
        rslt = modelComponents[["SEIR_model"]]$sample(samples, acceptFraction, 
                                                      batchSize)
        if (verbose) cat("Simulation complete\n")

        epsilon = rslt$result
        params = rslt$params

        cnames = paste("Beta_SE_", 
                       1:ncol(exposureModelInstance$X), sep = "")
        if (hasReinfection)
        {
            cnames = c(cnames, paste("Beta_RS_", 
                                     1:ncol(reinfectionModelInstance$X_prs), 
                                     sep = "")
            )
        }
        if (hasSpatial)
        {
            cnames = c(cnames, 
                       paste("rho_", 
                            1:length(distanceModelInstance$distanceList), 
                             sep = "")
            )
        }
        cnames = c(cnames, "gamma_EI", "gamma_IR")
        colnames(params) = cnames 
        
        modelResults[["param.samples"]] = params
        modelResults[["epsilon"]] = epsilon
        modelResults[["modelComponents"]] = list(
                     data_model = dataModelInstance,
                     exposure_model = exposureModelInstance,
                     reinfection_model = reinfectionModelInstance,
                     distance_model = distanceModelInstance,
                     transition_priors = transitionPriorsInstance,
                     initial_value_container = initialValueContainerInstance,
                     sampling_control = samplingControlInstance
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

plot.SpatialSEIRModel = function(modelObject)
{
    for (i in 1:ncol(modelObject$param.samples))
    {
        cat ("Press return to see next plot:")
        hist(modelObject$param.samples[,i], 
             main = colnames(modelObject$param.samples)[i],
             xlab = "Parameter Value",
             ylab = "MPD",
             freq=FALSE
             )
        if (i != ncol(modelObject$param.samples)) line <- readline()
    }


}

summary.SpatialSEIRModel = function(modelObject)
{
    nLoc = ncol(modelObject$modelComponents$data_model$Y)
    nTpt = nrow(modelObject$modelComponents$data_model$Y)

    hasSpatial = (modelObject$modelComponents$exposure_model$nLoc > 1) 
    hasReinfection = 
        (modelObject$modelComponents$reinfection_model$integerMode != 3) 

    qtiles = t(apply(modelObject$param.samples, 2, quantile, 
                   probs = c(0.025, 0.975)))

    means = apply(modelObject$param.samples, 2, mean)
    sds = apply(modelObject$param.samples, 2, sd)

    outMatrix = cbind(means, sds, qtiles)
    rownames(outMatrix) = colnames(modelObject$param.samples)
    colnames(outMatrix) = c("Mean", "SD", "95% LB", "95% UB")
    structure(list(modelComponents=modelObject$modelComponents,
       parameterEstimates=outMatrix,
       nLoc = nLoc,
       nTpt = nTpt,
       exposureParams = ncol(modelObject$modelComponents$exposure_model$X),
       reinfectionParams = ifelse(hasReinfection, 
            ncol(modelObject$modelComponents$reinfection_model$X_prs), 
            0),
       spatialParams = ifelse(hasSpatial, 
            length(modelObject$modelComponents$distance_model$distanceList),
            0)
       ), class = "summary.SpatialSEIRModel")
}


print.summary.SpatialSEIRModel = function(summaryObject)
{
    nl = "\n\n"
    cat(paste("Summary: SEIR Model", nl)) 
    cat(paste("Locations: ", summaryObject$nLoc, "\n",
              "Time Points: ", summaryObject$nTpt,"\n", sep = ""))
    cat(paste("Exposure Process Parameters: ", 
              summaryObject$exposureParams, "\n", sep = ""))
    cat(paste("Reinfection Model Parameters: ", 
             (summaryObject$reinfectionParams), "\n", sep = ""))
    cat(paste("Spatial Parameters: ", 
              (summaryObject$spatialParams), "\n", sep = ""))
    cat(nl)
    cat("Parameter Estimates:\n")
    print(round(summaryObject$parameterEstimates, 3))
    cat("\n") 
}
