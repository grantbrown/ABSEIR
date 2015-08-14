epidemic.simulations = function(modelObject, replicates=1, verbose = FALSE)
{
    if (class(modelObject) != "SpatialSEIRModel")
    {
        stop("modelObject must be of type SpatialSEIRModel")
    }

    modelCache = list()
    modelResult = list()
    tryCatch({
        dataModelInstance = modelObject$modelComponents$data_model;
        exposureModelInstance = modelObject$modelComponents$exposure_model;
        reinfectionModelInstance = 
            modelObject$modelComponents$reinfection_model;
        distanceModelInstance = modelObject$modelComponents$distance_model;
        transitionPriorsInstance = 
            modelObject$modelComponents$transition_priors;
        initialValueContainerInstance = 
            modelObject$modelComponents$initial_value_container; 
        samplingControlInstance = modelObject$modelComponents$sampling_control

        if (verbose) cat("...Building data model\n")
        modelCache[["dataModel"]] = new(dataModel, dataModelInstance$Y,
                                            dataModelInstance$type,
                                            dataModelInstance$compartment,
                                            dataModelInstance$phi,
                                            dataModelInstance$na_mask)


        if (verbose) cat("...Building distance model\n")
        modelCache[["distanceModel"]] = new(distanceModel)
        for (i in 1:length(distanceModelInstance$distanceList))
        {
            modelCache[["distanceModel"]]$addDistanceMatrix(
                distanceModelInstance$distanceList[[i]]
            )
        }
        modelCache[["distanceModel"]]$setPriorParameters(
            distanceModelInstance$priorAlpha,
            distanceModelInstance$priorBeta
        )

        if (verbose) cat("...Building exposure model\n")
        modelCache[["exposureModel"]] = new(
            exposureModel, 
            exposureModelInstance$X,
            exposureModelInstance$nTpt,
            exposureModelInstance$nLoc,
            exposureModelInstance$betaPriorMean,
            exposureModelInstance$betaPriorPrecision
        )
        if (!all(is.na(exposureModelInstance$offset)))
        {
            modelCache[["exposureModel"]]$offsets = (
                exposureModelInstance$offset
            )
        }

        if (verbose) cat("...Building initial value container\n")
        modelCache[["initialValueContainer"]] = new(initialValueContainer)
        modelCache[["initialValueContainer"]]$setInitialValues(
            initialValueContainerInstance$S0,
            initialValueContainerInstance$E0,
            initialValueContainerInstance$I0,
            initialValueContainerInstance$R0
        )

        if (verbose) cat("...Building reinfection model\n") 
        modelCache[["reinfectionModel"]] = new(
            reinfectionModel, 
            reinfectionModelInstance$integerMode
        )
        if (reinfectionModelInstance$integerMode != 3)
        {
            modelCache[["reinfectionModel"]]$buildReinfectionModel(
                reinfectionModelInstance$X_prs, 
                reinfectionModelInstance$priorMean, 
                reinfectionModelInstance$priorPrecision
        );
        }

        if (verbose) cat("...Building sampling control model\n") 
        modelCache[["samplingControl"]] = new (
            samplingControl, 
            c(samplingControlInstance$sim_width, samplingControlInstance$seed,
              samplingControlInstance$n_cores,samplingControlInstance$algorithm, 
              samplingControlInstance$batch_size,samplingControlInstance$epochs, 
              samplingControlInstance$max_batches),
            c(samplingControlInstance$acceptance_fraction, 
              samplingControlInstance$shrinkage)
        )

        if (verbose) cat("...building transition priors\n") 
        modelCache[["transitionPriors"]] = new(transitionPriors)
        modelCache[["transitionPriors"]]$setPriorsFromProbabilities(
            transitionPriorsInstance$p_ei,
            transitionPriorsInstance$p_ir,
            transitionPriorsInstance$p_ei_ess,
            transitionPriorsInstance$p_ir_ess
        )

        if (verbose) cat("Running epidemic simulations\n") 
       
        modelCache[["SEIRModel"]] = new( 
                    spatialSEIRModel, 
            modelCache[["dataModel"]],
            modelCache[["exposureModel"]],
            modelCache[["reinfectionModel"]],
            modelCache[["distanceModel"]],
            modelCache[["transitionPriors"]],
            modelCache[["initialValueContainer"]],
            modelCache[["samplingControl"]]
        )

        params = modelObject$param.samples
        params = params[rep(1:nrow(params), each = replicates),]
        modelResult[["simulatedResults"]] = 
            modelCache$SEIRModel$simulate(params)
        },
        warning=function(w){
            cat(paste("Warnings produced: ", w, sep = ""))
        },
        error=function(e){
            cat(paste("Errors produced: ", e, sep = ""))
        },
        finally={
            rm(modelCache)
        }
    );    
    names(modelResult$simulatedResults) = 
          paste("Simulation_", 1:length(modelResult$simulatedResults), 
                sep = "")
    return(modelResult$simulatedResults)
}


