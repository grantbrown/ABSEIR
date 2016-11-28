#' Update a \code{\link{SpatialSEIRModel}} object by drawing additional samples.
#'
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
    samples <- Ifelse("samples" %in% names(optparams), optparams$samples,  
                     nrow(modelObject$param.samples))
    sampling_control <- Ifelse("sampling_control" %in% names(optparams),
                              optparams$sampling_control, NA)
    verbose <- Ifelse("verbose" %in% names(optparams), 
                     optparams$verbose, FALSE)
    previous_epochs <- modelObject$completedEpochs  

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

        return(SpatialSEIRModel(dataModelInstance,
                                exposureModelInstance,
                                reinfectionModelInstance,
                                distanceModelInstance,
                                transitionPriorsInstance,
                                initialValueContainerInstance,
                                samplingControlInstance,
                                samples,
                                verbose,
                                is.updating = TRUE,
                                particles = modelObject$param.samples,
                                previous_eps = modelObject$current_eps,
                                previous_epochs = previous_epochs,
                                weights = modelObject$weights,
                                particle_eps = modelObject$epsilon))

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

