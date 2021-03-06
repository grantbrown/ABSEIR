#' Build a TransitionPriors object, which governs how individuals move from
#' the exposed to infectious and infectious to removed compartments. 
#' 
#' @param mode  The type of transition model to employ ("exponential" or "path_specific) 
#' @param params  Additional parameters, specific to the type of transition model. See
#' details section for additional information. 
#' 
#' @details 
#'  The TransitionPriors component of spatial SEIR(S) models captures the  
#'  process by which individuals move from the exposed to infectious compartment, 
#'  and from the infectious to removed compartment. This component thus governs the
#'  duration of the latent and infectious periods of the disease of interest, on the
#'  discrete timescale employed.
#' 
#'  Two different TransitionPriors configurations are currently offered: the exponential
#'  compartment membership model, and a path specific generalization. These may be
#'  specified manually using the TransitionPriors function, or with the specific
#'  ExponentialTransitionPriors and PathSpecificTransitionPriors functions. 
#' 
#'  The exponential version of this process requires four parameters: 
#'  \itemize{
#'      \item{p_ei}{The probability, at a given time point, that an exposed individual will 
#'                           become infectious} 
#'      \item{p_ir}{The probability, at a given time point, that an infectious individual 
#'                  be removed from the infectious population}
#'      \item{p_ei_ess}{An effective number of samples, corresponding to the confidence in
#'                      the chosen E to I transition probability}
#'      \item{p_ir_ess}{An effective number of samples, corresponding to the confidence in
#'                      the chosen I to R transition probability}
#'  }
#'  The path specific formulation requires fewer parameters, but more care is required in
#'  their specification. 
#'  \itemize{
#'   \item{Z1}{A probability density function for the time individuals spend in the latent state.} 
#'   \item{Z2}{A probability density function for the time individuals spend in the infectious state.} 
#'  }
#' 
#' @examples transitionPriors <- TransitionPriors("exponential", params=list(p_ei=
#'                                         1/5, p_ir=1/7, p_ei_ess=100, p_ir_ess=100))
#' 
#' @examples transitionPriors <- TransitionPriors("path_specific", params = list(
#'                                         Z1 = function(x){dunif(x, 2, 10)},
#'                                         Z2 = function(x){dunif(x, 7, 24)}))
#' @references "A path-specific SEIR model for use with general 
#'                latent and infectious time distributions." 2013. Porter, Aaron T, Oleson, Jacob J. Biometrics 69(1)
#' @export 
TransitionPriors = function(mode = c("exponential", "weibull", "path_specific"), params = list())
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
            return(structure(list(mode="exponential",
                       p_ei=params$p_ei,
                       p_ir=params$p_ir,
                       p_ei_ess=params$p_ei_ess,
                       p_ir_ess=params$p_ir_ess,
                       priorAlpha_gammaEI=NA,
                       priorBeta_gamma_EI=NA,
                       priorAlpha_gammaIR=NA,
                       priorBeta_gamma_IR=NA), class = "TransitionPriors"))
        }
    }
    else if (mode == "weibull")
    {
        requiredParams = c("latent_shape_prior_alpha", 
                           "latent_shape_prior_beta",
                           "latent_scale_prior_alpha",
                           "latent_scale_prior_beta",
                           "infectious_shape_prior_alpha",
                           "infectious_shape_prior_beta",
                           "infectious_scale_prior_alpha",
                           "infectious_scale_prior_beta")
        if (!all(sapply(requiredParams, function(x){x %in% names(params)})))
        {
            strParams = paste(requiredParams, collapse = ", ")
            stop(paste(
                       "The Weibull transition model ",
                       "requires the following parameters: ",
                       strParams, "\n", sep = ""))
        }

        minEIShape = qgamma(0.1, params$latent_shape_prior_alpha,
                                   params$latent_shape_prior_beta)
        maxEIScale = qgamma(0.9, params$latent_scale_prior_alpha,
                                   params$latent_scale_prior_beta)
        max_EI_idx = qweibull(1-1e-4, minEIShape, maxEIScale)

        minIRShape = qgamma(0.1, params$infectious_shape_prior_alpha,
                                   params$infectious_shape_prior_beta)
        maxIRScale = qgamma(0.9, params$infectious_scale_prior_alpha,
                                   params$infectious_scale_prior_beta)
        max_IR_idx = qweibull(1-1e-4, minIRShape, maxIRScale)

        if (max_EI_idx > 10000)
        {
            stop(paste("Your E to I transition prior allows individuals",
                          "to remain latent for up to ", max_EI_idx,  
                          ". Not recommended."))
        }
        else if (max_EI_idx > 100)
        {
            warning(paste("Your E to I transition prior allows individuals",
                          "to remain latent for up to ", max_EI_idx,  
                          ". This may be less efficient."))
        }
        
        if (max_IR_idx > 10000)
        {
            stop(paste("Your I to R transition prior allows individuals",
                          "to remain infectious for up to ", max_EI_idx,  
                          ". Not recommended."))
        }
        else if (max_IR_idx > 100)
        {
            warning(paste("Your I to R transition prior allows individuals",
                          "to remain infectious for up to ", max_EI_idx,  
                          ". This may be less efficient."))
        }
 

        return(structure(list(mode="weibull",
                   latent_shape_prior_alpha = params$latent_shape_prior_alpha,
                   latent_shape_prior_beta = params$latent_shape_prior_beta,
                   latent_scale_prior_alpha = params$latent_scale_prior_alpha,
                   latent_scale_prior_beta = params$latent_scale_prior_beta,
                   infectious_shape_prior_alpha = params$infectious_shape_prior_alpha,
                   infectious_shape_prior_beta = params$infectious_shape_prior_beta,
                   infectious_scale_prior_alpha = params$infectious_scale_prior_alpha,
                   infectious_scale_prior_beta = params$infectious_scale_prior_beta,
                   max_EI_idx = max_EI_idx,
                   max_IR_idx = max_IR_idx), 
                         class = "TransitionPriors"))
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

        Z1 = params$Z1
        Z2 = params$Z2
        # we can't assume that the functions passed in by the user are vectorized
        f1 = function(x){ sapply(x, Z1) }
        f2 = function(x){ sapply(x, Z2) }       


        inf.mean = integrate(function(x){x*f2(x)}, 0, Inf)$value
        hzd = function(x){
            sapply(x, Z2)/sapply(x, function(i){
                            1-integrate(f2, 0, i)$value
        })} 

        pdists = lapply(list(f1, f2), function(x){ 
            n = 100
            itrs = 0
            cmProb = 0
            # inefficient, but shouldn't be problematic
            while (itrs < 10 && cmProb < 1-truncation_prob)
            {
                indices = cbind(0:n,1:(n+1)) 
                probs = apply(indices, 1, function(x){integrate(f1, x[1], x[2])$value})
                cmProb = sum(probs)
                n = n*2
                itrs = itrs + 1
            }
            cprobs = cumsum(probs)
            sprobs = 1-c(0,cprobs[1:(length(cprobs)-1)])
            dhaz_probs = probs/sprobs
            chaz_probs = hzd(apply(indices,1,mean)) 
            max_idx = which(cprobs > 1-truncation_prob)[1]
            pdist = cbind(indices, probs, cprobs, sprobs, dhaz_probs,
                          chaz_probs)[1:max_idx,]
            colnames(pdist) = c("startTime", "endTime", "PDF",
                                "CDF", "Surv", "DiscHaz", "ContHaz")
            pdist
        })

        return(structure(list(mode="path_specific",
                              ei_pdist = pdists[[1]],
                              ir_pdist = pdists[[2]],
                              inf_mean = inf.mean), 
               class = "TransitionPriors"))
    }
    else
    {
        stop(paste("Transition mode: ", mode, " not recognized. Must be ",
             "exponential or path_specific\n", sep = ""))
    }
}

#' Build a path specific TransitionPriors object, which governs how individuals move from
#' the exposed to infectious and infectious to removed compartments. 
#' 
#' @param Z1 a probability density function for the length of time individuals spend in 
#'      the latent state
#' @param Z2 a probability density function for the length of time individuals spend in 
#'      the infectious state
#' @param truncation_prob ABSEIR precomputes the maximum number of discrete time points
#' during which a person may remain in the latent and infectious states. We therefore
#' require a probability cutoff. The trunction probability is the minimum that the
#' conditional probability of remaining in a disease state can become before individuals
#' are forced to transition. 
#' 
#' @details 
#'  The TransitionPriors component of spatial SEIR(S) models captures the  
#'  process by which individuals move from the exposed to infectious compartment, 
#'  and from the infectious to removed compartment. This component thus governs the
#'  duration of the latent and infectious periods of the disease of interest, on the
#'  discrete timescale employed.
#' 
#' @examples transitionPriors <- PathSpecificTransitionPriors(Z1 = function(x){dunif(x, 2, 10)},
#'                                         Z2 = function(x){dunif(x, 7, 24)})
#' @seealso \code{\link{TransitionPriors}}, \code{\link{ExponentialTransitionPriors}}
#' \code{\link{WeibullTransitionPriors}}
#' @export 
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
#' @seealso \code{\link{TransitionPriors}}, \code{\link{PathSpecificTransitionPriors}}, 
#' \code{\link{WeibullTransitionPriors}}
#' @export
ExponentialTransitionPriors = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    return(TransitionPriors("exponential", list(p_ei=p_ei,
                                                p_ir=p_ir,
                                                p_ei_ess=p_ei_ess,
                                                p_ir_ess=p_ir_ess)))
}


#' Build a path specific Weibull TransitionPriors object, with 
#' independent gamma priors for the shape/scale parameters of the 
#' latent and infectious time distributions. 
#' @param latent_shape_prior_alpha The alpha (shape) value for the 
#'    gamma hyperprior on the latent period shape parameter 
#' @param latent_shape_prior_beta The beta (rate) value for the gamma 
#'    hyperprior on the latent period shape parameter 
#' @param latent_scale_prior_alpha The alpha (shape) value for the gamma 
#'    hyperprior on the latent period scale parameter 
#' @param latent_scale_prior_beta The beta (rate) value for the gamma 
#'    hyperprior on the latent period scale parameter 
#' @param infectious_shape_prior_alpha The alpha (shape) value for the gamma 
#'    hyperprior on the infectious period shape parameter 
#' @param infectious_shape_prior_beta  The beta (rate) value for the gamma 
#'    hyperprior on the infectious period shape parameter 
#' @param infectious_scale_prior_alpha The alpha (shape) value for the gamma 
#'    hyperprior on the infectious period scale parameter 
#' @param infectious_scale_prior_beta The beta (rate) value for the gamma 
#'    hyperprior on the infectious period scale parameter 
#' @details 
#' WeibullTransitionPriors assumes that both the latent and infectious times
#' are distributed according to Weibull distributions. Let \eqn{\tau_{lat}} and
#' \eqn{\tau_{inf}} represent random variables corresponding to latent and infectious 
#' times. We then have: 
#' \eqn{\tau_{lat} \sim Weibull(\theta_1^{lat}, \theta_2^{lat})}
#' \eqn{\tau_{inf} \sim Weibull(\theta_1^{inf}, theta_2^{inf})}
#' For each of these four parameters, we assign a gamma prior, parameterized by 
#' shape and rate:
#' \eqn{f(\theta_{.}^{.}| \alpha, \beta^{-1}) = frac{1}{\Gamma(\alpha) \beta^\alpha}e^{-\frac{x}{\beta}}x^{\alpha-1}}
#'
#' These hyperparameters are specified by the eight arguments to WeibullTransitionPriors.
#'
#' 
#' @examples transitionPriors <- WeibullTransitionPriors(35,20,150,20,90,20,220,20)
#' @seealso \code{\link{TransitionPriors}}, \code{\link{PathSpecificTransitionPriors}}, 
#' \code{\link{ExponentialTransitionPriors}}
#' @export

WeibullTransitionPriors = function(latent_shape_prior_alpha, 
                                   latent_shape_prior_beta,
                                   latent_scale_prior_alpha,
                                   latent_scale_prior_beta,
                                   infectious_shape_prior_alpha,
                                   infectious_shape_prior_beta,
                                   infectious_scale_prior_alpha,
                                   infectious_scale_prior_beta)
{
    return(TransitionPriors("weibull", list(
                   latent_shape_prior_alpha = latent_shape_prior_alpha, 
                   latent_shape_prior_beta = latent_shape_prior_beta,
                   latent_scale_prior_alpha = latent_scale_prior_alpha,
                   latent_scale_prior_beta = latent_scale_prior_beta,
                   infectious_shape_prior_alpha = infectious_shape_prior_alpha,
                   infectious_shape_prior_beta = infectious_shape_prior_beta,
                   infectious_scale_prior_alpha = infectious_scale_prior_alpha,
                   infectious_scale_prior_beta = infectious_scale_prior_beta)))

}
