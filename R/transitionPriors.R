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

TransitionPriors = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    structure(list(mode="prob",
                   p_ei=p_ei,
                   p_ir=p_ir,
                   p_ei_ess=p_ei_ess,
                   p_ir_ess=p_ir_ess,
                   priorAlpha_gammaEI=NA,
                   priorBeta_gamma_EI=NA,
                   priorAlpha_gammaIR=NA,
                   priorBeta_gamma_IR=NA), class = "TransitionPriors")
}

