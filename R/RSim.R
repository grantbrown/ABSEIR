#' Simulate from a SEIR/SEIRS process using R (debug feature)
#' 
#' @param seed  A random seed
#' @param beta_SE  Chosen values of S to E transition coefficients
#' @param beta_RS  Chosen values of R to S transition coefficients
#' @param rho  Chosen values of distance parameters
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
#' @return a list
#' @details
#' This function is provided for testing and debugging purposes: it's not intended
#' for the end user
#'  
#' @examples \dontrun{results = Rsim(seed, 
#'                                    beta_SE,
#'                                    beta_RS, 
#'                                    rho,
#'                                    exposure_model,
#'                                    reinfection_model,
#'                                    distance_model,
#'                                    transition_priors,
#'                                    initial_value_container)}
#' @import Rcpp
#' @import methods
#' @useDynLib ABSEIR
#' @export
Rsim <- function(seed, 
                 beta_SE,
                 beta_RS, 
                 rho,
                 exposure_model,
                 reinfection_model,
                 distance_model,
                 transition_priors,
                 initial_value_container)
{ 
  nLoc <- exposure_model$nLoc
  nTpt <- exposure_model$nTpt
  X_SE <- exposure_model$X
  offs <- exposure_model$offset
  if (!all(is.na(offs))){
    stop("Temporal offsets not currently implemented")
  }

  hasSpatial <- ncol(distance_model$distanceList[[1]]) > 1
  nLags <- length(distance_model$laggedDistanceList[[1]])
  hasTSSpatial <- nLags > 0
  hasReinfection <- (reinfection_model$integerMode != 3)
  X_RS <- reinfection_model$X_prs
  
  if (transition_priors$mode != "exponential"){
    stop("Debugging currently only available for exponential transition modes")
  }
  gamma_EI <- -log(1-transition_priors$p_ei)
  gamma_IR <- -log(1-transition_priors$p_ir)
  p_EI <- transition_priors$p_ei
  p_IR <- transition_priors$p_ir
  
  set.seed(seed)

  E0 = initial_value_container$E0 
  I0 = initial_value_container$I0 
  R0 = initial_value_container$R0
  S0 = initial_value_container$S0
  N = E0 + I0 + R0 + S0
  
  DMlist = distance_model$distanceList
  lDMlist = distance_model$laggedDistanceList
    
  S=E=I=R=S_star=E_star=I_star=R_star=matrix(NA, nrow = nTpt, ncol = length(N))
  S[1,] = S0
  E[1,] = E0
  I[1,] = I0
  R[1,] = R0
  if (hasReinfection)
  {
    eta_RS <- exp(X_RS %*% beta_RS)
    p_RS <- 1-exp(-eta_RS)
  }
  eta_SE = matrix(exp(X_SE %*% beta_SE), nrow = nTpt, ncol = length(N))
  
  for (i in 1:nTpt)
  {
    intensity = (I[i,]/N * eta_SE[i,] +
                   apply(rho[1:length(DMlist)]*t(sapply(DMlist, function(x){
                     x %*% (I[i,]/N * eta_SE[i,])  
                   })), 2, sum))
    
    if (hasTSSpatial && i != 1)
    {
      lag.idx = 1
      for (j in (i-1):(max(i-length(lDMlist[[1]]), 1)))
      {
        intensity = intensity + rho[lag.idx+length(DMlist)]*(lDMlist[[j]][[lag.idx]] %*% (I[j,]/N*eta_SE[j,]))
        lag.idx = lag.idx + 1
      }
    }
    
    if (hasReinfection)
    {
      S_star[i,] = rbinom(n = rep(1, ncol(S_star)), 
                              size = R[i,],
                              prob = rep(p_RS[i], ncol(S_star)))
    }
    else
    {
      S_star[i,] = 0
    }
    
    p_SE = 1-exp(-intensity)
    E_star[i,] = rbinom(n = rep(1, ncol(E_star)), 
                        size = S[i,],
                        prob = p_SE)
    I_star[i,] = rbinom(n = rep(1, ncol(I_star)),
                        size = E[i,],
                        prob = p_EI)
    R_star[i,] = rbinom(n = rep(1, ncol(R_star)),
                        size = I[i,],
                        prob = p_IR)
    
    if (i != nTpt)
    {
      S[i+1,] = S[i,] + S_star[i,] - E_star[i,]
      E[i+1,] = E[i,] + E_star[i,] - I_star[i,]
      I[i+1,] = I[i,] + I_star[i,] - R_star[i,]
      R[i+1,] = R[i,] + R_star[i,] - S_star[i,]
    }
  }
  
  list(S=S,
       E=E,
       I=I,
       R=R,
       S0=S0,
       E0=E0,
       I0=I0,
       R0=R0,
       S_star=S_star,
       E_star=E_star,
       I_star=I_star,
       R_star=R_star)
} 
