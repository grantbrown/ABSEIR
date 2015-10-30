#' Create a SamplingControl object, which determines which ABC algorithm is 
#' to be used, and how it is configured. 
#' 
#' @param seed  an integer, giving the seed to be used when simulating epidemics
#' @param n_cores  an integer giving the number of CPU cores to employ
#' @param algorithm  a string, either equal to "BasicABC" for the simple
#' ABC rejection algorithm of Rubin (1980), or "Beaumont2009" for the
#' SMC approach of Beaumont et al. (2009). 
#' @param params optional algorithm configuration parameters, see: detail.  
#' @return an object of type \code{\link{SamplingControl}}
#' @details
#'  The basic ABC algorithm is useful in cases where good prior information is
#' available, and proposals from the prior distribution are therefore likely to 
#' fall with relative frequency into high posterior density regions. In cases
#' where the prior is diffuse with respect to the posterior, this can be very
#' inefficient. The SMC algorithm of Beaumont et al. 2009 randomly generates
#' parameters from a sequence of approximations to the posterior distribution, and
#' can greatly improve efficiency.\\
#' 
#' Additional parameters which may be passed to the algorithms:
#' \itemize{ 
#' \item{acceptance_fraction: }{For the BasicABC algorithm, this gives
#' the proportion of simulated epidemics to accept. The smaller the acceptance
#' fraction, the better the approximation to the posterior distribution, and the 
#' more computation time required.}
#' \item{batch_size: }{For the both algorithms, this determines the number of
#' epidemics to simulate in parallel, before returning to the main process to evaluate
#' them. \code{batch_size} must be greater than the number of samples requested 
#' in the \code{\link{SpatialSEIRModel}} function.}
#' \item{epochs: }{For the Beaumont2009 algorithm, \code{epochs} determines the maximum
#' number of iterations.}
#' \item{shrinkage: }{for the Beaumont2009 algorithm, \code{shrinkage} defines the multiplicative
#' constant by which the maximum distance between simulated and observed
#' epidemics is shrunk between each iteration.}
#' \item{max_batches: }{for the Beaumont2009 algorithm, \code{max_batches} determines
#' the maximum number of parallel batches to run before which a new set of 
#' parameters must be accepted. If an insufficient number of parameters are accepted
#' by the time the algorithm reaches \code{max_batches}, the program will terminate
#' under the assumption that the parameters have converged.}
#' \item{multivariate_perturbation}{A logical value indicating whether, for the
#' Beaumont2009 algorithm, parameter perturbations should be made from a
#' mulivariate normal distribuion rather than independent normals.}}
#' 
#' @examples samplingControl <- SamplingControl(123123, 2)
#' @export
SamplingControl = function(seed, n_cores, algorithm="Beaumont2009",
                           params=NA)                           
{
    alg = ifelse(algorithm == "BasicABC", 1,
          ifelse(algorithm == "Beaumont2009", 2, 
                 NA))
    if (is.na(alg))
    {
        stop(paste("Only the basic rejection algorithm and a modified version",
                   "of the SMC algorithm from Beaumont 2009 are currently",  
                   "supported.\n", sep = ""))
    }

    # Set default parameters where needed
    if (class(params) == "logical" && all(is.na(params))){
        if (algorithm == "Beaumont2009")
        { 
            params = list(acceptance_fraction = -1,
                 batch_size = 5000,
                 epochs = 100, 
                 shrinkage = 0.9,
                 max_batches = 20,
                 multivariate_perturbation = 0)
        }
        else
        {
            params = list(acceptance_fraction = 0.01,
                 batch_size = 10000,
                 epochs = 1,
                 shrinkage = 1,
                 max_batches = 1,
                 multivariate_perturbation = 0)
        }
    }
    else if (class(params) == "list")
    {
        if (algorithm == "Beaumont2009")
        {
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 5000
            }
            if (!("epochs" %in% names(params))) {
                params[["epochs"]] = 100
            }
            if (!("shrinkage" %in% names(params))) {
                params[["shrinkage"]] = 0.9
            }
            if (!("acceptance_fraction" %in% names(params))) {
                params[["acceptance_fraction"]] = 1
            }
            if (!("max_batches" %in% names(params))) {
                params[["max_batches"]] = 20
            }
            if (!("multivariate_perturbation" %in% names(params))){
                params[["multivariate_perturbation"]] = 0
            }
        }
        else if (algorithm == "BasicABC")
        {
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 10000
            }
            if (!("epochs" %in% names(params))) {
                params[["epochs"]] = 1
            }
            if (!("shrinkage" %in% names(params))) {
                params[["shrinkage"]] = 1
            }
            if (!("acceptance_fraction" %in% names(params))) {
                params[["acceptance_fraction"]] = 0.01
            }
            if (!("max_batches" %in% names(params))) {
                params[["max_batches"]] = 1
            }
            if (!("multivariate_perturbation" %in% names(params))){
                params[["multivariate_perturbation"]] = 0
            }
        }
        else
        {
            stop("params must be a list.")
        }
    }


    structure(list("sim_width" = 1,
                   "seed" = seed,
                   "n_cores" = n_cores,
                   "acceptance_fraction" = params$acceptance_fraction,
                   "batch_size" = params$batch_size,
                   "algorithm" = alg,
                   "epochs" = params$epochs,
                   "shrinkage" = params$shrinkage,
                   "max_batches" = params$max_batches,
                   "multivariate_perturbation" = params$multivariate_perturbation
                   ), class = "SamplingControl")
}


