#'  Compute empirically adjusted reproductive number curves for a PosteriorSimulation object
#'
#' @param  SimObject a PosteriorSimulation object, as created by the \code{\link{epidemic.simulations}} function.
#' @param cores  Optional argument - use multiple cores?
#'
#' @details  The main SpatialSEIRModel functon performs many simulations, but for the sake of
#'    memory efficiency and runtime does not return the simulated compartment values
#'    to the user. If simulated epidemics are desired, they may be quickly and
#'    easily generated using the \code{\link{epidemic.simulations}} function.
#'    Reproductive number estimation is performed using this function, as the
#'    calculations are somewhat computationally intensive and may not be required by
#'    all users.
#'
#' @examples \dontrun{r0 <- ComputeR0(epidemic.simulations(modelObject, replicates = 10,
#'                                                  verbose = TRUE))}
#'
#' @import parallel
#' @importFrom compiler cmpfun
#' @export
ComputeR0 <- function(SimObject, cores = 1)
{
  checkArgument("SimObject", mustHaveClass("PosteriorSimulation"))
  #cores <- SimObject$modelObject$modelComponents$sampling_control$n_cores
  cl <- makeCluster(cores)
  clusterExport(cl, "SimObject", envir = environment())
  setupR0 <- function(x){
      MO <<- SimObject$modelObject$modelComponents
      exposure_model <<- MO$exposure_model
      nTpt <<- exposure_model$nTpt
      nLoc <<- exposure_model$nLoc
      reinfection_model <<- MO$reinfection_model
      distance_model <<- MO$distance_model
      transition_priors <<- MO$transition_priors
      initial_value_container <<- MO$initial_value_container

      hasSpatial <<- ncol(distance_model$distanceList[[1]]) > 1
      nLags <<- length(distance_model$laggedDistanceList[[1]])
      hasTSSpatial <<- nLags > 0
      hasReinfection <<- (reinfection_model$integerMode != 3)
      X_RS <<- reinfection_model$X_prs
      X_SE <<- exposure_model$X
      nBetaRS <<- ifelse(hasReinfection, ncol(X_RS), 0)
      nBetaSE <<- ncol(X_SE)
      nRho <<- ifelse(hasSpatial, length(distance_model$distanceList), 0)
      nLRho <<- ifelse(hasTSSpatial, nLags, 0)

      DMlist <<- distance_model$distanceList
      lDMlist <<- distance_model$laggedDistanceList
  }
  parLapply(cl, 1:cores, setupR0)

  MO <- SimObject$modelObject$modelComponents
  exposure_model <- MO$exposure_model
  nTpt <- exposure_model$nTpt
  nLoc <- exposure_model$nLoc
  reinfection_model <- MO$reinfection_model
  distance_model <- MO$distance_model
  transition_priors <- MO$transition_priors
  initial_value_container <- MO$initial_value_container

  hasSpatial <- ncol(distance_model$distanceList[[1]]) > 1
  nLags <- length(distance_model$laggedDistanceList[[1]])
  hasTSSpatial <- nLags > 0
  hasReinfection <- (reinfection_model$integerMode != 3)
  X_RS <- reinfection_model$X_prs
  X_SE <- exposure_model$X
  nBetaRS <- ifelse(hasReinfection, ncol(X_RS), 0)
  nBetaSE <- ncol(X_SE)
  nRho <- ifelse(hasSpatial, length(distance_model$distanceList), 0)
  nLRho <- ifelse(hasTSSpatial, nLags, 0)

  DMlist <- distance_model$distanceList
  lDMlist <- distance_model$laggedDistanceList

  r0Func <- function(sim){
    paramvec <- SimObject$params[sim,]
    beta_SE <- paramvec[grepl("Beta_SE", names(paramvec), fixed = TRUE)]
    if (hasReinfection){
      beta_RS <- paramvec[grepl("Beta_RS", names(paramvec), fixed = TRUE)]
    }
    else{
      beta_RS <- c()
    }
    rho <- paramvec[grepl("rho", names(paramvec), fixed = TRUE)]

    if (transition_priors$mode == "exponential"){
        gamma_EI <- paramvec[grepl("gamma_EI", names(paramvec), fixed = TRUE)]
        gamma_IR <- paramvec[grepl("gamma_IR", names(paramvec), fixed = TRUE)]

        p_EI <- 1-exp(-gamma_EI)
        p_IR <- 1-exp(-gamma_IR)
    }
    else if (transition_priors$mode == "weibull"){
      gamma_EI <- NA
      gamma_IR <- NA

      EI_shape <- paramvec[names(paramvec) == "latent_shape"]
      EI_scale <- paramvec[names(paramvec) == "latent_scale"]

      IR_shape <- paramvec[names(paramvec) == "infectious_shape"]
      IR_scale <- paramvec[names(paramvec) == "infectious_scale"]

      n <- ceiling(qweibull(1-1e-4, shape = IR_shape, scale = IR_scale))
      indices <- cbind(0:n, 1:(n+1))
      p_IR_path <- c(apply(indices, 1, function(x){
        a <- pweibull(x[1], IR_shape, IR_scale)
        b <- pweibull(x[2], IR_shape, IR_scale)
        (b-a)/(1-a)
      }))
      p_IR_path <- p_IR_path[1:n]
    }
    else if (transition_priors$mode == "path-specific")
    {
      warning("Reproductive number estimation for general path-specific models unfinished")
      return(SimObject)
    }
    # Begin Calculation

    S0 = SimObject$simulationResults[[sim]]$S[1]
    E0 = SimObject$simulationResults[[sim]]$E[1]
    I0 = SimObject$simulationResults[[sim]]$I[1]
    R0 = SimObject$simulationResults[[sim]]$R[1]
    N = E0 + I0 + R0 + S0

    S <- SimObject$simulationResults[[sim]]$S
    E <- SimObject$simulationResults[[sim]]$E
    I <- SimObject$simulationResults[[sim]]$I
    R <- SimObject$simulationResults[[sim]]$R
    S_star <- SimObject$simulationResults[[sim]]$S_star
    E_star <- SimObject$simulationResults[[sim]]$E_star
    I_star <- SimObject$simulationResults[[sim]]$I_star
    R_star <- SimObject$simulationResults[[sim]]$R_star

    if (hasReinfection)
    {
      eta_RS <- exp(X_RS %*% beta_RS)
      p_RS <- 1-exp(-eta_RS)
    }
    eta_SE = matrix(exp(X_SE %*% beta_SE), nrow = nTpt, ncol = length(N))
    instantaneousExpectation <- matrix(NA, nrow = nTpt, ncol = nLoc)
    for (i in 1:nTpt)
    {
      # Non-spatial Component
      nonSpatialExpectation <- ifelse(I[i,] == 0, 0, S[i,]*(1-exp(-I[i,]/N * eta_SE[i,]))/I[i,])

      # Spatial Component
      k <- 1
      if (hasSpatial){
        spatialExpectation <- ifelse(I[i,] == 0, 0,
             apply(S[i,]*(1-exp(-t(rho[k]*(I[i,]/N * eta_SE[i,]) * t(DMlist[[k]])))),2,sum)/I[i,])
      }
      else{
        spatialExpectation <- rep(0, length(nonSpatialExpectation))
      }

      # Lagged component
      if (hasTSSpatial){
        laggedSpatialExpectation <- apply(sapply(1:length(lDMlist[[i]]), function(lag){
          if (lag + i > nTpt){
            return((ifelse(I[i,] == 0, 0,
                    apply(S[nTpt,]*(1-exp(-t(rho[nRho + lag]*(I[i,]/N * eta_SE[i,]) *
                                               t(lDMlist[[i]][[lag]])))),2,sum)/I[i,])))
          }
          (ifelse(I[i,] == 0, 0,
                  apply(S[i+lag,]*(1-exp(-t(rho[nRho + lag]*(I[i,]/N * eta_SE[i,]) *
                                              t(lDMlist[[i]][[lag]])))),2,sum)/I[i,]))
        }),1,sum)
      }
      else{
        laggedSpatialExpectation <- rep(0, length(nonSpatialExpectation))
      }

      instantaneousExpectation[i,] <- nonSpatialExpectation + spatialExpectation + laggedSpatialExpectation
    }
    r_EA <- instantaneousExpectation
    if (transition_priors$mode == "exponential"){
      p_IR <- 1-exp(-gamma_EI)
      for (i in 1:nrow(r_EA))
      {
        pI <- 1
        idx <- i
        while (pI > 1e-8)
        {
          pI <- pI*(1-p_IR)
          idx <- ifelse(idx == nrow(r_EA), idx, idx+1)
          r_EA[i,] <- r_EA[i,] + instantaneousExpectation[idx,]*pI
        }
      }
    }
    else if (transition_priors$mode == "weibull"){
      for (i in 1:nrow(r_EA))
      {
        pI <- 1
        idx <- i
        surv <- 1
        while (pI > 1e-8 && surv <= length(p_IR_path))
        {
          pI <- pI*(1-p_IR_path[surv])
          surv <- surv + 1
          idx <- ifelse(idx == nrow(r_EA), idx, idx+1)
          r_EA[i,] <- r_EA[i,] + instantaneousExpectation[idx,]*pI
        }
      }
    }
    else if (transition_priors$mode == "path_specific"){
      # To-do
    }
    colnames(r_EA) <- paste("location_", 1:ncol(r_EA), sep = "")
    r_EA
  }

  r0FuncC <- cmpfun(r0Func)
  repNums <- parLapply(cl,
            1:length(SimObject$simulationResults),
            r0FuncC)

  stopCluster(cl)
  print("Done with R0sim")
  for (sim in 1:length(repNums)){
    SimObject$simulationResults[[sim]][["R_EA"]] <- repNums[[sim]]
  }
  return(SimObject)

}

