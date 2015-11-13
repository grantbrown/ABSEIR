#' perform and return epidemic simulations based on a fitted model object
#' 
#' @param modelObject a SpatialSEIRModel object, as created by the \code{\link{SpatialSEIRModel}}
#' function. 
#' @param replicates the number of replicate simulations to perform for each sample
#' contained in \code{modelObject}
#' @param verbose a logical value, indicating whether verbose output should be 
#' provided. 
#' @param returnCompartments  a logical value, indicating whether simulated compartments
#' should be returned along with epsilon values
#' 
#' @details 
#'    The main SpatialSEIRModel functon performs many simulations, but for the sake of 
#'    memory efficiency and runtime does not return the simulated compartment values
#'    to the user. If simulated epidemics are desired, they may be quickly and
#'    easily generated using this function. 
#' 
#' @examples \dontrun{simulate_values <- epidemicSimulations(modelObject, replicates = 10, 
#'                                                  returnCompartments = TRUE, 
#'                                                  verbose = TRUE)} 
#' 
#' @export
epidemic.simulations = function(modelObject, replicates=1, returnCompartments = TRUE,verbose = FALSE)
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
                                            dataModelInstance$cumulative,
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
              samplingControlInstance$max_batches, 
              samplingControlInstance$multivariate_perturbation),
            c(samplingControlInstance$acceptance_fraction, 
              samplingControlInstance$shrinkage)
        )

        if (verbose) cat("...building transition priors\n") 
        modelCache[["transitionPriors"]] = new(transitionPriors,
                                               transitionPriorsInstance$mode)
        transitionMode = transitionPriorsInstance$mode
        if (transitionMode == "exponential")
        {
            modelCache[["transitionPriors"]]$setPriorsFromProbabilities(
                transitionPriorsInstance$p_ei,
                transitionPriorsInstance$p_ir,
                transitionPriorsInstance$p_ei_ess,
                transitionPriorsInstance$p_ir_ess
            )
        }
        else if (transitionMode == "weibull")
        {
            modelCache[["transitionPriors"]]$setPriorsForWeibull(
                              c(transitionPriorsInstance$latent_shape_prior_alpha,
                                transitionPriorsInstance$latent_shape_prior_beta,
                                transitionPriorsInstance$latent_scale_prior_alpha,
                                transitionPriorsInstance$latent_scale_prior_beta),
                              c(transitionPriorsInstance$infectious_shape_prior_alpha,
                                transitionPriorsInstance$infectious_shape_prior_beta,
                                transitionPriorsInstance$infectious_scale_prior_alpha,
                                transitionPriorsInstance$infectious_scale_prior_beta),
                                transitionPriorsInstance$max_EI_idx,
                                transitionPriorsInstance$max_IR_idx)
        }
        else
        {
            modelCache[["transitionPriors"]]$setPathSpecificPriors(
                                            transitionPriorsInstance$ei_pdist,
                                            transitionPriorsInstance$ir_pdist,
                                            transitionPriorsInstance$inf_mean)
        }

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
        if (returnCompartments)
        {
            modelResult[["simulatedResults"]] = 
                modelCache$SEIRModel$simulate(params)
        } 
        else
        {
            modelResult[["simulatedResults"]] = 
                modelCache$SEIRModel$evaluate(params)
        }
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
    if (returnCompartments)
    {
        names(modelResult$simulatedResults) = 
              paste("Simulation_", 1:length(modelResult$simulatedResults), 
                    sep = "")
    }
    return(structure(list(modelObject = modelObject, 
                          simulationResults=modelResult$simulatedResults),
                     class = "PosteriorSimulation"))
}


