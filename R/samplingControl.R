#' Create a SamplingControl object, which determines which ABC algorithm is 
#' to be used, and how it is configured. 
#' 
#' @param seed  an integer, giving the seed to be used when simulating epidemics
#' @param n_cores  an integer giving the number of CPU cores to employ
#' @param algorithm  a string, either equal to "BasicABC" for the simple
#' ABC rejection algorithm of Rubin (1980), "Beaumont2009" for the
#' SMC approach of Beaumont et al. (2009), or "DelMoral2006"  for the adaptive
#' SMC approach of Del Moral (2012). 
#' @param  params optional algorithm configuration parameters, see: detail.  
#' @return an object of type \code{\link{SamplingControl}}
#' @details
#'  The basic ABC algorithm is useful in cases where good prior information is
#' available, and proposals from the prior distribution are therefore likely to 
#' fall with relative frequency into high posterior density regions. In cases
#' where the prior is diffuse with respect to the posterior, this can be very
#' inefficient. The SMC algorithms of Del Moral 2012 and
#' Beaumont et al. 2009 randomly generate
#' parameters from a sequence of approximations to the posterior distribution, and
#' can greatly improve efficiency.\\
#' 
#' Additional parameters which may be passed to the algorithms:
#' \itemize{ 
#' \item{acceptance_fraction: }{For the BasicABC algorithm, this gives
#' the proportion of simulated epidemics to accept. The smaller the acceptance
#' fraction, the better the approximation to the posterior distribution, and the 
#' more computation time required.}
#' \item{target_eps:}{For all algorithms, this determines an epsilon value at which 
#' the program will terminate, declaring convergence.}
#' \item{batch_size: }{For all algorithms, this determines the number of
#' epidemics to simulate in parallel, before returning to the main process to evaluate
#' them. \code{batch_size} must be greater than the number of samples requested 
#' in the \code{\link{SpatialSEIRModel}} function.}
#' \item{init_batch_size: }{For all algorithms, this optionally determines
#' a distinct batch size for the first iteration.}
#' \item{epochs: }{For the Beaumont2009 and DelMoral2012 algorithms, 
#' \code{epochs} determines the maximum
#' number of iterations.}
#' \item{shrinkage: }{for the Beaumont2009 algorithm, \code{shrinkage} defines the multiplicative
#' constant by which the maximum distance between simulated and observed
#' epidemics is shrunk between each iteration. For the DelMoral2012 algorithm, this
#' parameter determines the quality index, \eqn{\alpha}{alpha} between zero and one.}
#' \item{max_batches: }{for the Beaumont2009 and DelMoral2012 algorithms, \code{max_batches} determines
#' the maximum number of parallel batches to run before which a new set of 
#' parameters must be accepted. If an insufficient number of parameters are accepted
#' by the time the algorithm reaches \code{max_batches}, the program will terminate
#' under the assumption that the parameters have converged.}
#' \item{multivariate_perturbation}{A logical value indicating whether, for the
#' Beaumont2009 algorithm, parameter perturbations should be made from a
#' mulivariate normal distribuion rather than independent normals.}
#' \item{m}{For the DelMoral2012 algorithm, an integer determining the number of 
#' simulated epidemics to use for each set of basis parameters (parameterized
#' as in the 2012 paper.)}
#' \item{particles}{For the 'simulate' algorithm, a raw matrix of parameter 
#' values must be provided, following the format created by the other
#' algorithms. This functonality is used primarily for debugging purposes; most
#' users should perform such simulations using the 
#' \code{\link{epidemic.simulations}} function instead.}
#' \item{replicates}{For the 'simulate' algorithm, a number of replicate
#' simulations to be performed per particle.}}
#' 
#' 
#' @examples samplingControl <- SamplingControl(123123, 2)
#' @export
SamplingControl = function(seed, n_cores, algorithm="Beaumont2009",
                           params=NA)                           
{
    alg = ifelse(algorithm == "BasicABC", 1,
          ifelse(algorithm == "Beaumont2009", 2, 
          ifelse(algorithm == "DelMoral2012", 3, 
          ifelse(algorithm == "simulate", 4, 
                 NA))))
    if (is.na(alg) || alg == "DelMoral2012" || alg == "BasicABC")
    {
        stop(paste("Algorithm", algorithm, "not supported. Choice must be one of:",
                   "Beaumont2009"))
    }

    # Set default parameters where needed
    if (class(params) == "logical" && all(is.na(params))){
        if (algorithm == "Beaumont2009")
        { 
            params = list(acceptance_fraction = -1,
                 target_eps = 0,
                 batch_size = 5000,
                 init_batch_size = 5000,
                 epochs = 100, 
                 shrinkage = 0.9,
                 max_batches = 20,
                 multivariate_perturbation = 0,
                 m=1,
                 particles=-1,
                 replicates=-1)
        }
        else if (algorithm == "DelMoral2012")
        {
             params = list(acceptance_fraction = -1,
                 target_eps = 0,
                 batch_size = 5000,
                 init_batch_size = 5000,
                 epochs = 100, 
                 shrinkage = 0.9,
                 max_batches = 20,
                 multivariate_perturbation = 0,
                 m=5,
                 particles=-1,
                 replicates=-1)           
        }
        else if (algorithm == "simulate")
        {
            stop("the 'simulate' algorithm requires pre-specified parameters.")
        }
        else
        {
            params = list(acceptance_fraction = 0.01,
                 target_eps = 0,
                 batch_size = 10000,
                 init_batch_size = 10000,
                 epochs = 1,
                 shrinkage = 1,
                 max_batches = 1,
                 multivariate_perturbation = 0,
                 m=1,
                 particles=-1,
                 replicates=-1)
        }
    }
    else if (class(params) == "list")
    {   
        # Some parameters were provided, set defaults for any which were not
        if (algorithm == "Beaumont2009")
        {
            if (!("target_eps" %in% names(params))){
                params[["target_eps"]] = 0
            }
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 5000
            }
            if (!("init_batch_size" %in% names(params))) {
                params[["init_batch_size"]] = params[["batch_size"]]
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
            if (!("m" %in% names(params))){
                params[["m"]] = 1
            }
            if (!("particles" %in% names(params))){
                params[["particles"]] = -1
            }
            if (!("replicates" %in% names(params))){
                params[["replicates"]] = -1
            }

        }
        else if (algorithm == "DelMoral2012")
        {
            if (!("target_eps" %in% names(params))){
                params[["target_eps"]] = 0
            }
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 5000
            }
            if (!("init_batch_size" %in% names(params))) {
                params[["init_batch_size"]] = params[["batch_size"]]
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
            if (!("m" %in% names(params))){
                params[["m"]] = 5
            }
            if (!("particles" %in% names(params))){
                params[["particles"]] = -1
            }
            if (!("replicates" %in% names(params))){
                params[["replicates"]] = -1
            }
        }
        else if (algorithm == "simulate")
        {
            if (!("target_eps" %in% names(params))){
                params[["target_eps"]] = 0
            }
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 5000
            }
            if (!("init_batch_size" %in% names(params))) {
                params[["init_batch_size"]] = params[["batch_size"]]
            }
            if (!("epochs" %in% names(params))) {
                params[["epochs"]] = 0
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
            if (!("m" %in% names(params))){
                params[["m"]] = 1
            }
            if (!("particles" %in% names(params))){
                stop("The 'particles' argument is required for the 'simulate' algorithm")
            }
            if (!("replicates" %in% names(params))){
                params[["replicates"]] = 1
            }
        }
        else if (algorithm == "BasicABC")
        {
            if (!("target_eps" %in% names(params))){
                params[["target_eps"]] = 0
            }
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 10000
            }
            if (!("init_batch_size" %in% names(params))) {
                params[["init_batch_size"]] = params[["batch_size"]]
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
            if (!("m" %in% names(params))){
                params[["m"]] = 1
            }
            if (!("particles" %in% names(params))){
                params[["particles"]] = -1
            }
            if (!("replicates" %in% names(params))){
                params[["replicates"]] = -1
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
                   "target_eps" = params$target_eps,
                   "batch_size" = params$batch_size,
                   "init_batch_size" = params$init_batch_size,
                   "algorithm" = alg,
                   "epochs" = params$epochs,
                   "shrinkage" = params$shrinkage,
                   "max_batches" = params$max_batches,
                   "multivariate_perturbation" = params$multivariate_perturbation,
                   "m"=params$m,
                   "particles"=params$particles,
                   "replicates"=params$replicates
                   ), class = "SamplingControl")
}


