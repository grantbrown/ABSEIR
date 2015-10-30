#' compute approximate Bayes Factor in favor of one spatial SEIR model over another
#' 
#' @param model1 A model in favor of which one would like to compute a Bayes Factor 
#' function. 
#' @param model2 The model to which \code{model1} is to be compared 
#' @param p1 The prior probability of \code{model1} 
#' @param p2 The prior probability of \code{model2} 
#' @param n_samples The desired number of accepted simulation values on which the 
#' Bayes Factor calculation is to be based
#' @param batch_size The number of epidemics to simulate in parallel before 
#' assessing the number of accepted samples
#' @param max_itrs The maximum number of parallel batches to execute before giving up 
#' 
#' @details A Bayes Factor is a measure of the posterior evidence in favor
#' of one model compared to another. In the ABC setting, we may compute  
#' approximate Bayes Factors of comparably converged models by assessing
#' the parameter acceptance rate at a new iteration. 
#' 
#'
#' @examples \dontrun{compareModels(model1, model2)}
#'                                                
#' @export 
compareModels = function(model1, model2, p1 = NA, p2 = NA, n_samples = 1000,
                         batch_size = 10000, max_itrs = 1000)
{
    if (class(model1) != "SpatialSEIRModel" || 
        class(model2) != "SpatialSEIRModel")
    {
       stop("compareModels may only be used to compare SpatialSEIRModel objects.") 
    }
    else if (any(dim(model1$modelComponents$data_model$Y) !=
                 dim(model1$modelComponents$data_model$Y)) ||
             any(model1$modelComponents$data_model$Y != 
                 model1$modelComponents$data_model$Y))
    {
       stop(paste(
            "Only models of the same data should be compared in this fashion, ",
            "or at all really...", sep = "")) 
    }
    else
    {
        earlyTerm1 = (model1$completedEpochs != 
                      model1$modelComponents$sampling_control$epochs)
        earlyTerm2 = (model2$completedEpochs != 
                      model2$modelComponents$sampling_control$epochs)
        match1 = (model1$modelComponents$sampling_control$max_batches == 
                  model1$modelComponents$sampling_control$max_batches)
        match2 = (model1$modelComponents$sampling_control$batch_size ==
                  model2$modelComponents$sampling_control$batch_size)
        test1 = !(earlyTerm1 && earlyTerm2 && match1 && match2)
        test2 = abs(model1$current_eps - model2$current_eps) > 1e-8  
        if (test1 && test2)
        {
            warning(paste(
                    "Models to be compared should either terminate at the ",
                    "same epsilon value, or represent the same terminating ",
                    "acceptance rate. Comparison may be unreliable."))
        }     
    }

    if (all(is.na(p1)) && all(is.na(p2)))
    {
        p1 = 0.5
        p2 = 0.5
        print("Assuming equal prior probabilities.")
    }
    else if (all(is.na(p1) != all(is.na(p2))))
    {
        stop("If some prior probabilities are supplied, all must be supplied.")
    }



    w1 = model1$weights
    w2 = model1$weights
    e1 = model1$current_eps
    e2 = model2$current_eps
    e.compare = max(c(e1, e2))
    
    drawSamples = function()
    {
        mr1 = model1
        mr2 = model2
        s1 = model1$param.samples[sample(1:nrow(model1$param.samples), 
                                         prob = w1, replace = TRUE,
                                         size = batch_size),]
        s2 = model2$param.samples[sample(1:nrow(model1$param.samples), 
                                         prob = w1, replace = TRUE,
                                         size = batch_size),]
        mr1$param.samples = s1
        mr2$param.samples = s2
        e1 = epidemic.simulations(mr1, returnCompartments=FALSE
                                  )$simulationResults$result
        e2 = epidemic.simulations(mr2, returnCompartments=FALSE
                                  )$simulationResults$result
        c(sum(e1 < e.compare), sum(e2 < e.compare))
    }

    itrs = 0
    accepted = c(0,0)
    while (itrs < max_itrs && sum(accepted) < n_samples)
    {
        itrs = itrs + 1
        accepted = accepted + drawSamples() 
    }
    if (itrs == 1000)
    {
        warning("n_samples not reached before max_itrs")
    }

    accepted = accepted/(itrs*batch_size)
    out = accepted[1]/accepted[2]
    names(out) = "Approximate Bayes Factor: model1 vs. model2"
    out
}
