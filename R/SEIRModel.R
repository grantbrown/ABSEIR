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

    modelComponents = list()
    if (verbose){cat(paste("Initializing Model Components", sep = ""))}
    result = tryCatch({
        if (verbose) cat(paste("...Building data model", sep = ""))
        modelComponents[["dataModel"]] = new(dataModel, dataModelInstance$Y,
                                            dataModelInstance$type,
                                            dataModelInstance$compartment,
                                            dataModelInstance$phi)

        if (verbose) cat(paste("...Building distance model", sep = ""))
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

        if (verbose) cat(paste("...Building exposure model", sep = ""))
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

        if (verbose) cat(paste("...Building initial value container", sep = ""))
        modelComponents[["initialValueContainer"]] = new(initialValueContainer)
        modelComponents[["initialValueContainer"]]$setInitialValues(
            initialValueContainerInstance$S0,
            initialValueContainerInstance$E0,
            initialValueContainerInstance$I0,
            initialValueContainerInstance$R0
        )

        if (verbose) cat(paste("...Building reinfection model", sep = "")) 
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

        if (verbose) cat(paste("...Building sampling control model", sep = "")) 
        modelComponents[["samplingControl"]] = new (
            samplingControl, 
            samplingControlInstance$sim_width,
            samplingControlInstance$seed,
            samplingControlInstance$n_cores
        )

        if (verbose) cat(paste("...Building transition priors", sep = "")) 
        modelComponents[["transitionPriors"]] = new(transitionPriors)
        modelComponents[["transitionPriors"]]$setPriorsFromProbabilities(
            transitionPriorsInstance$p_ei,
            transitionPriorsInstance$p_ir,
            transitionPriorsInstance$p_ei_ess,
            transitionPriorsInstance$p_ir_ess
        )

        if (verbose) cat(paste("...Preparing model object", sep = ""))
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
        if (verbose) cat(paste("Running main simulation", sep = ""))
        modelComponents[["SEIR_model"]]$sample(samples, acceptFraction, batchSize)
        if (verbose) cat(paste("Simulation complete", sep = ""))
    },
    warning=function(w)
    {
        cat(paste("Warnings generated:", sep = ""))
        print(w)
    },
    error=function(e)
    {
        cat(paste("Errors generated:", sep = ""))
        print(e)
        if (exists("modelComponents"))
        {
            rm(modelComponents)
        }
        stop("Aborting")
    },
    finally = { 
        if (exists("modelComponents"))
        {
            rm(modelComponents)
        }
    })
    return(structure(result, class = "SpatialSEIRModel"))
}

