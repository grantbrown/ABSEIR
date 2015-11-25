#' compute approximate Bayes Factor in favor of one spatial SEIR model over another
#' 
#' @param modelList A list of models to compare with approximate Bayes factors.  
#' function. 
#' @param priors The prior probabilities of each model in \code{modelList} 
#' @param n_samples The desired number of accepted simulation values on which the 
#' Bayes Factor calculation is to be based
#' @param batch_size The number of epidemics to simulate in parallel before 
#' assessing the number of accepted samples
#' @param max_itrs The maximum number of parallel batches to execute before giving up 
#' @param epsilon The cutoff value used to determine whether simulated epidemics 
#' are accepted or rejected. If left blank, the mean of the two smallest terminating
#' epsilon values models under comparison is used. If these are dramatically different,
#' this approach may produce misleading results. 
#' @param verbose A logical value, indicating whether progress information should
#' be displayed. 
#' 
#' @details A Bayes Factor is a measure of the posterior evidence in favor
#' of one model compared to another. In the ABC setting, we may compute  
#' approximate Bayes Factors of comparably converged models by assessing
#' the parameter acceptance rate at a new iteration. 
#' 
#'
#' @examples \dontrun{compareModels(list(model1, model2))}
#'                                                
#' @export 
compareModels = function(modelList, priors=NA, n_samples = 1000,
                         batch_size = 10000, max_itrs = 1000,
                         epsilon=NA, verbose=FALSE)
{
    correctClasses = sapply(modelList, function(x){class(x)  == 
                            "SpatialSEIRModel"})

    correctDim = sapply(modelList, function(x){
                        all(dim(x$modelComponents$data_model$Y) == 
                            dim(modelList[[1]]$modelComponents$data_model$Y)) &&
                        all(x$modelComponents$data_model$Y == 
                            modelList[[1]]$modelComponents$data_model$Y)})
    if (any(!correctClasses))
    {
       stop("compareModels may only be used to compare SpatialSEIRModel objects.") 
    }
    else if (any(!correctDim))
    {
       stop(paste(
            "Only models of the same data should be compared in this fashion, ",
            "or at all really...", sep = "")) 
    }
    else
    {
        earlyTerm = sapply(modelList, function(x){
                           x$completedEpochs != 
                                x$modelComponents$sampling_control$epochs})
        match = sapply(modelList, function(x){
                       x$modelComponents$sampling_control$max_batches ==
               modelList[[1]]$modelComponents$sampling_control$max_batches
                                })

        test1 = !(all(earlyTerm) && all(match))
        test2 = max(abs(sapply(modelList, function(x){
                x$current_eps - modelList[[1]]$current_eps
        }))) > 1e-8
        if (test1 && test2)
        {
            warning(paste(
                    "Models to be compared should either terminate at the ",
                    "same epsilon value, or represent the same terminating ",
                    "acceptance rate. Comparison may be unreliable."))
        }     
    }

    if (all(is.na(priors)))
    {
        priors = rep(1/length(modelList), length(modelList))
        print("Assuming equal prior probabilities.")
    }
    else if (length(priors) != length(modelList) || sum(priors) != 1)
    {
        stop("If some prior probabilities are supplied, all must be supplied, and must sum to 1.")
    }

    weightList = lapply(modelList, function(x){x$weights})
    epsVec = sapply(modelList, function(x){x$current_eps})

    if (is.na(epsilon))
    {
        e.compare = mean(epsVec[order(epsVec)][1:2])
    }
    else
    {
        e.compare = epsilon
    }
    
    drawSamples = function()
    {
        mr = modelList
        s = lapply(modelList, function(x){
                   x$param.samples[sample(1:nrow(x$param.samples),
                                          prob = x$weights, replace = TRUE,
                                          size = batch_size),]
                                })
        for (i in 1:length(s))
        {
            mr[[i]]$param.samples = s[[i]]
        }

        esim = lapply(1:length(mr), function(x){
                      if (verbose)
                      {
                          cat(paste("  Evaluating model ", x, "\n", sep = ""))
                      }
                      epidemic.simulations(mr[[x]], 
                        returnCompartments=FALSE)$simulationResults$result
                                })
        sapply(esim, function(x){
               sum(x < e.compare)})
    }

    itrs = 0
    accepted = rep(0, length(modelList))
    while (itrs < max_itrs && sum(accepted) < n_samples)
    {
        itrs = itrs + 1
        if (verbose)
        {
            cat(paste("Iteration ", itrs, "\n", sep = ""))
        }
        accepted = accepted + drawSamples() 
        if (verbose)
        {
            cat(paste(" ", sum(accepted), " of ", n_samples, " samples obtained.\n", sep = ""))
        }
    }
    if (itrs == max_itrs)
    {
        warning("n_samples not reached before max_itrs")
    }

    BF = accepted/(itrs*batch_size)
    BF = BF*priors
    return(matrix(BF, nrow = length(BF), ncol = length(BF)) /
    matrix(BF, nrow = length(BF), ncol = length(BF), byrow = TRUE))
}

