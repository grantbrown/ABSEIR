test_that("Spatial models with no contact structure produce correct sims",{

## Lagged contact test. 
    seed = 123412
    trueBeta = c(-1.7, -1.5)
    trueRho = c(0.2, 0.1, 0.02, 0.01, 0.01, 0.01)*0
    gamma_EI = 1/5
    gamma_IR = 1/7
    N = c(1000000, 1000000, 1000000, 1000000)
    E0 = c(100, 100, 100,100)
    I0 = c(10, 10, 10,10)
    R0 = c(0,0,0,0)
    S0 = N - E0 - I0 - R0

    dummy_I_star <- matrix(0, nrow = 100, ncol = length(S0))


    design_matrix = cbind(1,            # Intercept 
                          c(rep(0, 50), # Intervention
                            1:50/50))

    X1 <- cbind(1,            # Intercept 
                c(rep(0, 50), # Intervention
                  1:50/50))

    X2 <- cbind(1,            # Intercept 
                c(rep(0, nrow(X1)))) # Intervention

    X3 <- cbind(1,            # Intercept 
                c(rep(0, 10), # Intervention
                  1:90/50))

    X4 <- cbind(1,            # Intercept 
                c(rep(0, 90), # Intervention
                  1:10))

    design_matrix <- rbind(X1,X2,X3,X4)



#[rep(1:nrow(design_matrix), ncol(DMlist[[1]])),]

    DMlist = list(
      matrix(c(0,1,1,0,
               1,0,0,0,
               1,0,0,0,
               0,0,0,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,0,0,
               0,0,1,0,
               0,1,0,0,
               0,0,0,0), nrow = 4, byrow = TRUE)
    )
    lDMlist = lapply(1:nrow(dummy_I_star), function(x){
      list(
        matrix(c(0,0,0,0,
                 0,0,0,0,
                 0,0,0,0,
                 0,0,0,0), nrow = 4, byrow = TRUE),
        matrix(c(0,0,0,0,
                 0,0,0,0,
                 0,0,0,0,
                 0,0,0,0), nrow = 4, byrow = TRUE),
        matrix(c(0,0,0,0,
                 0,0,0,0,
                 0,0,0,0,
                 0,0,0,0), nrow = 4, byrow = TRUE),
        matrix(c(0,0,0,0,
                 0,0,0,0,
                 0,0,0,0,
                 0,0,0,0), nrow = 4, byrow = TRUE)
      )
    })

    lDMlist[[50]] = list(
      matrix(c(0,0,1,1,
               0,0,1,1,
               1,1,0,1,
               1,1,1,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,1,1,
               0,0,1,1,
               1,1,0,1,
               1,1,1,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,1,1,
               0,0,1,1,
               1,1,0,1,
               1,1,1,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,1,1,
               0,0,1,1,
               1,1,0,1,
               1,1,1,0), nrow = 4, byrow = TRUE)
    )


    data_model = DataModel(dummy_I_star,
                           type = "identity",
                           compartment="I_star",
                           cumulative=FALSE)

    exposure_model = ExposureModel(design_matrix,
                                     nTpt = nrow(dummy_I_star),
                                     nLoc = ncol(dummy_I_star),
                                     betaPriorPrecision = 0.5,
                                     betaPriorMean = 0)

    reinfection_model = ReinfectionModel("SEIR")


    distance_model = TDistanceModel(distanceList = DMlist,
                                    laggedDistanceList = lDMlist,
                                    priorAlpha = 1, 
                                    priorBeta = 10)


    initial_value_container = InitialValueContainer(S0=S0,
                                                    E0=E0,
                                                    I0=I0,
                                                    R0=R0)

    transition_priors = ExponentialTransitionPriors(p_ei = 1-exp(-1/5), 
                                                    p_ir= 1-exp(-1/7),
                                                    p_ei_ess = 100,
                                                    p_ir_ess = 100)

    partMat <- matrix(c(trueBeta, trueRho, gamma_EI, gamma_IR), nrow = 1)
    colnames(partMat) <- c(paste("Beta_SE_", 1:length(trueBeta), sep = ""), 
                           paste("rho_", 1:length(trueRho), sep = ""),
                           "gamma_EI", "gamma_IR")
      
    sampling_control = SamplingControl(seed = 123123, 
                                       n_cores = 1,
                                       algorithm="simulate",
                                       list(particles = partMat,
                                            replicates = 1000, 
                                            batch_size = 1000))

    runtime = system.time(result <- SpatialSEIRModel(data_model,
                                                       exposure_model,
                                                       reinfection_model,
                                                       distance_model,
                                                       transition_priors,
                                                       initial_value_container,
                                                       sampling_control,
                                                       samples = 5,
                                                       verbose = FALSE))

    runtime2 <- system.time(dat.m <- lapply(1:1000, function(x){
      Rsim(seed=123115+x, 
           beta_SE = trueBeta,
           beta_RS = c(), 
           rho = trueRho,
           exposure_model=exposure_model,
           reinfection_model=reinfection_model,
           distance_model=distance_model,
           transition_priors=transition_priors,
           initial_value_container=initial_value_container)
    }))

    sims <- result



    I_star1_R <- sapply(dat.m, function(x){x$I_star[,1]})
    I_star2_R <- sapply(dat.m, function(x){x$I_star[,2]})
    I_star3_R <- sapply(dat.m, function(x){x$I_star[,3]})
    I_star4_R <- sapply(dat.m, function(x){x$I_star[,4]})

    I_star1 <- sapply(sims$simulationResults, function(x){x$I_star[,1]})
    I_star2 <- sapply(sims$simulationResults, function(x){x$I_star[,2]})
    I_star3 <- sapply(sims$simulationResults, function(x){x$I_star[,3]})
    I_star4 <- sapply(sims$simulationResults, function(x){x$I_star[,4]})
    I1mean <- apply(I_star1, 1, median)
    I2mean <- apply(I_star2, 1, median)
    I3mean <- apply(I_star3, 1, median)
    I4mean <- apply(I_star4, 1, median)
    I1mean_R <- apply(I_star1_R, 1, median)
    I2mean_R <- apply(I_star2_R, 1, median)
    I3mean_R <- apply(I_star3_R, 1, median)
    I4mean_R <- apply(I_star4_R, 1, median)

    I1LB <- apply(I_star1, 1, quantile, probs = c(0.1))
    I2LB <- apply(I_star2, 1, quantile, probs = c(0.1))
    I3LB <- apply(I_star3, 1, quantile, probs = c(0.1))
    I4LB <- apply(I_star4, 1, quantile, probs = c(0.1))

    I1UB <- apply(I_star1, 1, quantile, probs = c(0.9))
    I2UB <- apply(I_star2, 1, quantile, probs = c(0.9))
    I3UB <- apply(I_star3, 1, quantile, probs = c(0.9))
    I4UB <- apply(I_star4, 1, quantile, probs = c(0.9))


    I1LB_R <- apply(I_star1_R, 1, quantile, probs = c(0.1))
    I2LB_R <- apply(I_star2_R, 1, quantile, probs = c(0.1))
    I3LB_R <- apply(I_star3_R, 1, quantile, probs = c(0.1))
    I4LB_R <- apply(I_star4_R, 1, quantile, probs = c(0.1))

    I1UB_R <- apply(I_star1_R, 1, quantile, probs = c(0.9))
    I2UB_R <- apply(I_star2_R, 1, quantile, probs = c(0.9))
    I3UB_R <- apply(I_star3_R, 1, quantile, probs = c(0.9))
    I4UB_R <- apply(I_star4_R, 1, quantile, probs = c(0.9))

    ckfc <- function(a,b){
      var(a) == var(b) ||
        cor(a, b) >= 0.95 &&
        mean(abs(a-b))/mean(b) < 0.1
    }


    passed <- TRUE
    if (!(ckfc(I1LB, I1LB_R) &&
      ckfc(I2LB, I2LB_R) &&
      ckfc(I3LB, I3LB_R) &&
      ckfc(I4LB, I4LB_R) &&
      ckfc(I1UB, I1UB_R) &&
      ckfc(I2UB, I2UB_R) &&
      ckfc(I3UB, I3UB_R) &&
      ckfc(I4UB, I4UB_R) &&
      ckfc(I1mean, I1mean_R) &&
      ckfc(I2mean, I2mean_R) &&
      ckfc(I3mean, I3mean_R) &&
      ckfc(I4mean, I4mean_R))){
        passed <- FALSE
    }
    expect(passed, "R and C++ simulations do not agree.")
})


