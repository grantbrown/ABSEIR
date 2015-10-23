#' Build a TransitionPriors object, which governs how individuals move from
#' the exposed to infectious and infectious to removed compartments. 
#' 
#' @param p_ei the estimated probability that, at a given time point, an exposed 
#' individual becomes infectious. 
#' @param p_ir the estimated probability that, at a given time point, an infectious 
#' individual is removed from the population.
#' @param p_ei_ess the effective sample size, or corresponding number of observations
#' of the E to I process, associated with \code{p_ei}
#' @param p_ir_ess the effective sample size, or corresponding number of observations
#' of the I to R process, associated with \code{p_ir}
#' 
#' @details 
#'  Many other transition models are possible, and more are planned for this package. 
#'  In particular, we would like to be able to support path specific transition
#'  trajectories for individuals in these compartments. 
#' 
#'  In the meantime, we use the parameterization from Brown et al. 2015, which 
#'  employs a constant transition probability over discrete time, corresponding to
#'  a geometric compartment membership time. For additional discussion, see Brown 
#'  et al. 2015. 
#' 
#' @examples transitionPriors <- TransitionPriors(1/5,1/7, 100, 100)
TransitionPriors = function(mode = c("exponential", "path_specific"), params = list())
{
    if (class(params) != "list")
    {
        stop("The params argument must be a list")
    }

    if (mode == "exponential")
    {
        if (!("p_ei" %in% names(params)) || 
            !("p_ir" %in% names(params)) || 
            !("p_ei_ess" %in% names(params)) || 
            !("p_ir_ess" %in% names(params)))
        {
            stop(paste("Exponential transition priors require four parameters: ", 
                 "p_ei, p_ir, p_ei_ess, and p_ir_ess", sep = ""))
        }
        else
        {
            return(ExponentialTransitionPriors(params$p_ei, params$p_ir,
                                               params$p_ei_ess, 
                                               params$p_ir_ess))
        }
    }
    else if (mode == "path_specific")
    {
        if (!("Z1" %in% names(params)) || 
            !("Z2" %in% names(params)))
        {
            stop(paste("Path specific transition models require a distribution ",
                       "for the latent period, Z1, and for the infectious period, Z2\n",
                       sep = ""))
        }
        truncation_prob = ifelse("truncation_prob" %in% names(params), 
                                 params$truncation_prob, 1e-6)
        if ()
        Z1 = params$Z1
        Z2 = params$Z2
        # we can't assume that the functions passed in by the user are vectorized
        f1 = function(x){ sapply(x, Z1) }
        f2 = function(x){ sapply(x, Z2) }       
    }
    else
    {
        stop(paste("Transition mode: ", mode, " not recognized. Must be ",
             "exponential or path_specific\n", sep = ""))
    }
}

PathSpecificTransitionPriors = function(Z1 = NA, Z2 = NA, truncation_prob=1e-6)
{
    if (length(Z1) != 1 || length(Z2) != 1 || 
        class(Z1) != "function" || class(Z2) != "function")
    {
       stop("Z1 and Z2 must be functions describing the latent and infectious time PDFs.") 
    }
    else
    {
        return(TransitionPriors("path_specific", list("Z1" = Z1,
                                                      "Z2" = Z2,
                                                      "truncation_prob" = 
                                                      truncation_prob)))
    }
}


#' Build an exponential TransitionPriors object, which governs how individuals move from
#' the exposed to infectious and infectious to removed compartments. 
#' 
#' @param p_ei the estimated probability that, at a given time point, an exposed 
#' individual becomes infectious. 
#' @param p_ir the estimated probability that, at a given time point, an infectious 
#' individual is removed from the population.
#' @param p_ei_ess the effective sample size, or corresponding number of observations
#' of the E to I process, associated with \code{p_ei}
#' @param p_ir_ess the effective sample size, or corresponding number of observations
#' of the I to R process, associated with \code{p_ir}
#' 
#' @details 
#'  ExponentialTransitionPriors provides the parameterization from Brown et al. 2015, which 
#'  employs a constant transition probability over discrete time, corresponding to
#'  a geometric compartment membership time. For additional discussion, see Brown 
#'  et al. 2015. 
#' 
#' @examples transitionPriors <- ExponentialTransitionPriors(1/5,1/7, 100, 100)
ExponentialTransitionPriors = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    structure(list(mode="exponential",
                   p_ei=p_ei,
                   p_ir=p_ir,
                   p_ei_ess=p_ei_ess,
                   p_ir_ess=p_ir_ess,
                   priorAlpha_gammaEI=NA,
                   priorBeta_gamma_EI=NA,
                   priorAlpha_gammaIR=NA,
                   priorBeta_gamma_IR=NA), class = "TransitionPriors")
}

