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
        rslt = modelComponents[["SEIR_model"]]$sample(samples, acceptFraction, batchSize)
        if (verbose) cat("Simulation complete\n")

        epsilon = rslt$result
        params = rslt$params

        # add names here 
        
        modelResults[["results"]] = rslt
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
    return(structure(modelResults[["results"]], class = "SpatialSEIRModel"))
}

