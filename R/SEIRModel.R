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
#'        specified by the sampling_control argument. Ignored for the 
#'        debug-oriented 'simulate' algorithm.
#' @param verbose print diagnostic information on the progress of the fitting algorithm.
#'        Available output levels are 0, 1, 2, and 3, in ascending order
#'        of detail. Level 0 output will
#'        print almost no progress/diagnostic information to the log. Level 1 
#'        Will provide iteration updates only. Leve 2 provides additional 
#'        chain setup and diagnostic information. Level 3 prints calculation
#'        diagnostic information. 
#' @param \dots Additional arguments, used internally. 
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
                          verbose=FALSE,
                          ...)
{
    checkArgument("samples", mustHaveClass(c("integer",
                                                      "numeric")),
                             mustBeLen(1),
                             mustBeInRange(lower = 1)
                          ) 
    # Todo: convert following checks to new infrastructure

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

    if (sampling_control$algorithm == 4)
    {
        # We don't actually want to fit a model
        sampling_control2 <- SamplingControl(seed = sampling_control$seed, 
                                             n_cores = sampling_control$n_cores)
        sampling_control2$epochs <- 0
        sampling_control2$batch_size <- 1
        dummy_model <- SpatialSEIRModel(data_model,
                                        exposure_model,
                                        reinfection_model,
                                        distance_model,
                                        transition_priors,
                                        initial_value_container,
                                        sampling_control2,
                                        samples=1,
                                        verbose=FALSE)
        dummy_model$modelComponents$sampling_control <- sampling_control
        dummy_model$param.samples <- sampling_control$particles
        return(epidemic.simulations(dummy_model, 
                                    sampling_control$replicates,
                                    verbose=verbose))
    }
    # Check if we're in update mode
    optionalParamNames <- c("particles", "is.updating", "previous_eps",
                            "previous_epochs", "weights", "particle_eps")
    optparams <- list(...)
    is.updating <- Ifelse("is.updating" %in% names(optparams), 
                         optparams$is.updating,  
                         FALSE)
    particles <- Ifelse("particles" %in% names(optparams),
                              optparams$particles, NA)
    previous_eps <- Ifelse("previous_eps" %in% names(optparams), 
                     optparams$previous_eps, NA)
    previous_epochs <- Ifelse("previous_epochs" %in% names(optparams),
                              optparams$previous_epochs, NA)
    previous_weights <- Ifelse("weights" %in% names(optparams), 
                               optparams$weights, NA)
    start_result <- Ifelse("particle_eps" %in% names(optparams), 
                           optparams$particle_eps, NA)
    if (is.updating && (any(is.na(previous_eps)) || any(is.na(particles)) 
                        || any(is.na(previous_epochs)) ||
                        any(is.na(previous_weights)) || 
                        any(is.na(start_result))))
    {
        stop("In update mode, particles, weights, epochs, eps vector, and previous epsilon required")
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
                                            c(data_model$phi,
                                              data_model$report_fraction,
                                              data_model$report_fraction_ess),
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
            if (exposure_model$nTpt != length(distance_model$laggedDistanceList))
            {
                stop("Lagged distance model and exposure model imply different number of time points.")
            }

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
        modelComponents[["initialValueContainer"]] = new(initialValueContainer, initial_value_container$type)
        modelComponents[["initialValueContainer"]]$setInitialValues(
            initial_value_container$S0,
            initial_value_container$E0,
            initial_value_container$I0,
            initial_value_container$R0,
            initial_value_container$max_S0,
            initial_value_container$max_E0,
            initial_value_container$max_I0,
            initial_value_container$max_R0
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
              sampling_control$batch_size,
              sampling_control$init_batch_size,
              sampling_control$epochs, 
              sampling_control$max_batches, 
              sampling_control$multivariate_perturbation, 
              sampling_control$m),
            c(sampling_control$acceptance_fraction, sampling_control$shrinkage,
              sampling_control$lpow,sampling_control$target_eps
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
        if (verbose && is.updating) cat("Initializing existing parameters.\n")
        if (is.updating)
        {
            modelComponents[["SEIR_model"]]$setParameters(particles, 
                                                          previous_weights, 
                                                          start_result,
                                                          previous_eps)
        }
        if (verbose) cat("Running main simulation\n")
        rslt = modelComponents[["SEIR_model"]]$sample(samples, 
                                                      sampling_control$keep_compartments*1, 
                                                      verbose*1)

        if (verbose) cat("Simulation complete\n")

        epsilon = rslt$result
        params = rslt$params
        weights = rslt$weights
        current_eps = rslt$currentEps
        completed_epochs = Ifelse(is.updating, 
                                  rslt$completedEpochs + previous_epochs, 
                                  rslt$completedEpochs)
                                  
        

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
        if (data_model$type == "fractional")
        {
            cnames <- c(cnames, "report_fraction")
        }
        cnames <- c(cnames, 
                    paste0("S0_", 1:length(initial_value_container$S0)),
                    paste0("E0_", 1:length(initial_value_container$E0)),
                    paste0("I0_", 1:length(initial_value_container$I0)),
                    paste0("R0_", 1:length(initial_value_container$R0)))


        if (length(cnames) == ncol(params)){
            colnames(params) = cnames 
        } else{
            print(cnames)
            print(dim(params))

            warning("Parameters are of unexpected dimension! Failed to name")
        }
        
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
        if (sampling_control$keep_compartments > 0){
            modelResults[["simulationResults"]] = rslt$simulationResults
        }
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


