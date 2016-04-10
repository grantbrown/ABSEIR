#' Fit a spatial/non-spatial SEIR/SEIRS model based on the provided model components.
#' 
#' @param data_model A data model object, describing the link between
#'  the observed data and the unobserved epidemic counts. Valid data models
#'  are created using the \code{\link{DataModel}} function. 
#' @param exposure_model An exposure model object, which describes the
#'  spatial and temporal variability of the exposure/infection process. Valid 
#'  exposure models are created using the \code{\link{ExposureModel}} function. 
#' @param reinfection_model A reinfection model object, which describes
#'  whether or not individuals are able to return from the Removed category
#'  to the Susceptible population. Valid reinfection models are created using the
#'  \code{\link{ReinfectionModel}} function. 
#' @param distance_model A distance model object which describes the underlying
#'  contact network, in addition to prior parameters which constrain the contact
#'  process. Valid distance models are created using \code{\link{DistanceModel}} 
#' @param transition_priors An object containing information about the E to I and
#'  I to R transition prior parameters. These are created using the 
#'  \code{\link{TransitionPriors}} function, and involve either transition
#'  probabilities and corresponding effective sample sizes (amount of prior 
#'  information), or manually specified probability distributions. .  
#' @param initial_value_container An object specifying the initial state
#'  of the epidemic for each spatial location, created by the
#'  \code{\link{InitialValueContainer}} function. 
#' @param sampling_control An object specifying information about the sampling
#'  algorithm. In particular, the sampling_control argument should specify the
#'  number of CPU cores to employ, and the random seed to use. 
#'  Sampling control objects are created by the 
#'  \code{\link{SamplingControl}} function. 
#' @param samples the number of samples to approximate from the posterior
#'        distribution, i.e. the number of particles to simulate. The number
#'        of particles should be considerably smaller than the batch size
#'        specified by the sampling_control argument. 
#' @param verbose print diagnostic information on the progress of the fitting algorithm
#' @return an object of type \code{\link{SpatialSEIRModel}} 
#' @details
#' Use the supplied model components to build and fit a corresponding model.   
#' This function is used to fit all of the models in the spatial SEIRS model class. 
#'    Numerous ABC algorithms have been developed, but as of now ABSEIR provides 
#'    just two. The first algorithm is the basic rejection algorithm of Rubin 
#'    1980. While this approach performs well when good prior information is available, 
#'    it can be extremely inefficient when prior distributions are diffuse with 
#'    respect to the posterior. To address this shortcoming, we have implemented
#'    the Sequential Monte-Carlo approach proposed by Beaumont 2009, 2010. We may
#'    provide additional algorithms in the future, in particular that of Del Moral
#'    et al. 2012. 
#'  
#' @examples \dontrun{results = SpatialSEIRModel(data_model, exposure_model,
#'                                                 reinfection_model, distance_model,
#'                                                 transition_priors, initial_value_container,
#'                                                 sampling_control, 50, TRUE)}
#' @seealso \code{\link{DataModel}}, \code{\link{ExposureModel}}, 
#' \code{\link{ReinfectionModel}}, \code{\link{DistanceModel}}, 
#' \code{\link{TransitionPriors}}, \code{\link{InitialValueContainer}}, 
#' \code{\link{SamplingControl}}, \code{\link{summary.SpatialSEIRModel}}, 
#' \code{\link{plot.SpatialSEIRModel}}, \code{\link{compareModels}}, 
#' \code{\link{epidemic.simulations}}, 
#' @import Rcpp
#' @import methods
#' @useDynLib ABSEIR
#' @export
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
    transitionMode= transition_priors$mode
    result = tryCatch({
        if (verbose) cat("...Building data model\n")
        modelComponents[["dataModel"]] = new(dataModel, data_model$Y,
                                            data_model$type,
                                            data_model$compartment,
                                            data_model$cumulative,
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
        nLags <- length(distance_model$laggedDistanceList[[1]]) 
        modelComponents[["distanceModel"]]$setupTemporalDistanceMatrices(
                                    exposure_model$nTpt
                                )
        if (nLags > 0)
        {
            for (i in 1:length(distance_model$laggedDistanceList)) 
            {
                for (j in 1:nLags)
                { 
                    modelComponents[["distanceModel"]]$addTDistanceMatrix(i,
                                distance_model$laggedDistanceList[[i]][[j]]
                    ) 
                }
            }
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
              sampling_control$max_batches, 
              sampling_control$multivariate_perturbation),
            c(sampling_control$acceptance_fraction, sampling_control$shrinkage,
              sampling_control$target_eps
              )
        )

        if (verbose) cat("...Building transition priors\n") 
        modelComponents[["transitionPriors"]] = new(transitionPriors, 
                                                    transition_priors$mode)
        if (transitionMode == "exponential")
        {
            modelComponents[["transitionPriors"]]$setPriorsFromProbabilities(
                transition_priors$p_ei,
                transition_priors$p_ir,
                transition_priors$p_ei_ess,
                transition_priors$p_ir_ess
            )
        }
        else if (transitionMode == "weibull")
        {
            modelComponents[["transitionPriors"]]$setPriorsForWeibull(
                              c(transition_priors$latent_shape_prior_alpha,
                                transition_priors$latent_shape_prior_beta,
                                transition_priors$latent_scale_prior_alpha,
                                transition_priors$latent_scale_prior_beta),
                              c(transition_priors$infectious_shape_prior_alpha,
                                transition_priors$infectious_shape_prior_beta,
                                transition_priors$infectious_scale_prior_alpha,
                                transition_priors$infectious_scale_prior_beta),
                                transition_priors$max_EI_idx,
                                transition_priors$max_IR_idx)
        }
        else
        {
            modelComponents[["transitionPriors"]]$setPathSpecificPriors(
                                                    transition_priors$ei_pdist,
                                                    transition_priors$ir_pdist,
                                                    transition_priors$inf_mean)
        }

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
        weights = rslt$weights
        current_eps = rslt$currentEps
        completed_epochs = rslt$completedEpochs

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
                            1:distance_model$len, 
                             sep = "")
            )
        }
        if (transitionMode == "exponential")
        { 
            cnames = c(cnames, "gamma_EI", "gamma_IR")
        }
        else if (transitionMode == "weibull")
        {
            cnames = c(cnames, 
                       "latent_shape", 
                       "latent_scale",
                       "infectious_shape",
                       "infectious_scale"
                       )

        }
        colnames(params) = cnames 
        
        modelResults[["param.samples"]] = params
        modelResults[["epsilon"]] = epsilon
        modelResults[["weights"]] = weights
        modelResults[["current_eps"]] = current_eps
        modelResults[["modelComponents"]] = list(
                     data_model = data_model,
                     exposure_model = exposure_model,
                     reinfection_model = reinfection_model,
                     distance_model = distance_model,
                     transition_priors = transition_priors,
                     initial_value_container = initial_value_container,
                     sampling_control = sampling_control
        ) 
        modelResults[["completedEpochs"]] = completed_epochs
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


#' Update a \code{\link{SpatialSEIRModel}} object by drawing additional samples.
#' @param object a \code{\link{SpatialSEIRModel}} object
#' @param \dots Additional arguments include:
#' \itemize{
#'  \item{\code{sampling_control}: }{a \code{\link{SamplingControl}} object}
#'  \item{\code{samples}: }{indicates the number of parameter values to sample 
#' from the posterior}
#'  \item{\code{verbose}: }{indicates whether verbose output should be 
#' provided.}
#' }
#' @export
update.SpatialSEIRModel = function(object, ...)
{
    modelObject = object
    if (class(modelObject) != "SpatialSEIRModel")
    {
        stop("object must be of type SpatialSEIRModel")
    }

    optionalParamNames <- c("samples", "sampling_control", "verbose")
    optparams <- list(...)
    samples = Ifelse("samples" %in% names(optparams), optparams$samples,  
                     nrow(modelObject$param.samples))
    sampling_control = Ifelse("sampling_control" %in% names(optparams),
                              optparams$sampling_control, NA)
    verbose = Ifelse("verbose" %in% names(optparams), 
                     optparams$verbose, FALSE)

    if (!(all(is.na(sampling_control))) && class(sampling_control) 
        != "SamplingControl")
    {
        stop("The sampling_control argument must be of type SamplingControl.")
    }
    else if ((class(sampling_control) == "SamplingControl") && 
             (sampling_control$algorithm != 
             modelObject$modelComponents$sampling_control$algorithm))
    {
        stop(paste("The algorithm specified by sampling_control must match ", 
             "the algorithm used in the original model.\n", sep = "")) 
    } 

    modelCache = list()
    modelResults = list()
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
        if (all(is.na(sampling_control))){
            samplingControlInstance = modelObject$modelComponents$sampling_control
        }
        else
        {
            samplingControlInstance = sampling_control 
        } 

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
        nLags <- length(distanceModelInstance$laggedDistanceList[[1]]) 
        modelCache[["distanceModel"]]$setupTemporalDistanceMatrices(
                    exposureModelInstance$nTpt
                ) 
        if (nLags > 0)
        {
            for (i in 1:length(distanceModelInstance$laggedDistanceList)) 
            {
                for (j in 1:nLags)
                {
                    modelCache[["distanceModel"]]$addTDistanceMatrix(i,
                                distanceModelInstance$laggedDistanceList[[i]][[j]]
                    ) 
                }
            }
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
              samplingControlInstance$shrinkage, 
              samplingControlInstance$target_eps)
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
                transitionPriorsInstance$p_ir_ess)
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


        if (verbose) cat("Running additional epidemic simulations\n") 
       
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

        rslt = modelCache[["SEIRModel"]]$update(samples, 
                                                modelObject$param.samples,
                                                modelObject$epsilon, 
                                                modelObject$weights,
                                                modelObject$current_eps,
                                                verbose*1)
        if (verbose) cat("Simulation complete\n")

        epsilon = rslt$result
        params = rslt$params
        weights = rslt$weights
        current_eps = rslt$currentEps
        completed_epochs = rslt$completedEpochs

        cnames = paste("Beta_SE_", 
                       1:ncol(exposureModelInstance$X), sep = "")
        hasSpatial = (ncol(dataModelInstance$Y) > 1) 
        hasReinfection = (reinfectionModelInstance$integerMode != 3) 

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
                            1:(distanceModelInstance$len), 
                             sep = "")
            )
        }

        if (transitionMode == "exponential")
        {
            cnames = c(cnames, "gamma_EI", "gamma_IR")
        }
        else if (transitionMode == "weibull")
        {
            cnames = c(cnames, 
                       "latent_shape", 
                       "latent_scale",
                       "infectious_shape",
                       "infectious_scale"
                       )

        }

        colnames(params) = cnames 
        
        modelResults[["param.samples"]] = params
        modelResults[["epsilon"]] = epsilon
        modelResults[["weights"]] = weights
        modelResults[["current_eps"]] = current_eps
        modelResults[["modelComponents"]] = list(
                     data_model = dataModelInstance,
                     exposure_model = exposureModelInstance,
                     reinfection_model = reinfectionModelInstance,
                     distance_model = distanceModelInstance,
                     transition_priors = transitionPriorsInstance,
                     initial_value_container = initialValueContainerInstance,
                     sampling_control = samplingControlInstance)
        modelResults[["completedEpochs"]] = completed_epochs
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
                rm(modelCache)
            }
        })
        return(structure(modelResults, class = "SpatialSEIRModel"))   
}

#' Plot a graphical summary of the marginal posterior distribution of SEIR model parameters
#' @param x a \code{\link{SpatialSEIRModel}} object
#' @param \dots additional arguments to be passed to plotting functions. 
#' @examples \dontrun{plot(modelObject)}
#' @export
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
#' @export
summary.SpatialSEIRModel = function(object, ...)
{
    nLoc = ncol(object$modelComponents$data_model$Y)
    nTpt = nrow(object$modelComponents$data_model$Y)

    transitionMode = object$modelComponents$transition_priors$mode
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
       reinfectionParams = Ifelse(hasReinfection, 
            ncol(object$modelComponents$reinfection_model$X_prs), 
            0),
       spatialParams = Ifelse(hasSpatial, 
            length(object$modelComponents$distance_model$distanceList),
            0),
       transitionParams = Ifelse(transitionMode == "exponential", 2,
                          Ifelse(transitionMode == "weibull", 4, 0))
       ), class = "summary.SpatialSEIRModel")
}

#' Print a \code{summary.SpatialSEIRModel} object
#' @param x a \code{summary.SpatialSEIRModel} object
#' @param \dots not used
#' @export
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
    cat(paste("Transition Parameters: ", 
              (x$transitionParams), "\n", sep = ""))

    cat(nl)
    cat("Parameter Estimates:\n")
    print(round(x$parameterEstimates, 3))
    cat("\n") 
}
