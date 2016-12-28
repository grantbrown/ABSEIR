test_that("Spatial models with lagged contact produce correct sims",{

    seed = 123412
    trueBeta = c(-0.5, -2)
    trueRho = c(0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25)
    gamma_EI = 1/5
    gamma_IR = 1/7
    N = c(1000000, 1000000, 1000000, 1000000)
    E0 = c(100,0,0,0)
    I0 = c(20,0,0,0)
    R0 = c(0,0,0,0)
    S0 = N - E0 - I0 - R0

    dummy_I_star <- matrix(0, nrow = 100, ncol = length(S0))

    X1 <- cbind(1,            # Intercept 
                c(rep(0, 20), # Intervention
                  1:80/80))

    X2 <- cbind(1,            # Intercept 
                c(rep(0, nrow(X1)))) # Intervention

    design_matrix <- rbind(X2,X2,X2,X2)

    DMlist = list(
      matrix(c(0,0,0,0,
               0,0,0,0,
               0,0,0,0,
               1,0,0,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,0,0,
               0,0,0,0,
               0,0,0,1,
               0,0,0,0), nrow = 4, byrow = TRUE),
      matrix(c(0,0,0,0,
               0,0,1,0,
               0,0,0,0,
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



    lDMlist[[10]] = list(
      matrix(c(0,0,0,0,
               1,0,0,0,
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
               1,0,0,0,
               0,0,0,0), nrow = 4, byrow = TRUE)
    )



    lDMlist[[30]] = list(
      matrix(c(0,0,0,0,
               1,0,0,0,
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
                                       n_cores = 10,
                                       algorithm="simulate",
                                       list(particles = partMat,
                                            replicates = 1000,
                                            batch_size = 1))


    runtime = system.time(result <- SpatialSEIRModel(data_model,
                                                       exposure_model,
                                                       reinfection_model,
                                                       distance_model,
                                                       transition_priors,
                                                       initial_value_container,
                                                       sampling_control,
                                                       samples = 1000,
                                                       verbose = FALSE))

    Rcomp <- lapply(1:1000, function(x){
      Rsim(seed=123115+x, 
           params = list(beta_SE = trueBeta,
                         beta_RS = c(), 
                         rho = trueRho,
                         gamma_EI = gamma_EI,
                         gamma_IR = gamma_IR),
           exposure_model=exposure_model,
           reinfection_model=reinfection_model,
           distance_model=distance_model,
           transition_priors=transition_priors,
           initial_value_container=initial_value_container)
    })

    CppSim <- result
    compareComp <- function(cmpName, Rcomp, CppSim){
      Comp_1_R <- sapply(Rcomp, function(x){x[[cmpName]][,1]})
      Comp_2_R <- sapply(Rcomp, function(x){x[[cmpName]][,2]})
      Comp_3_R <- sapply(Rcomp, function(x){x[[cmpName]][,3]})
      Comp_4_R <- sapply(Rcomp, function(x){x[[cmpName]][,4]})
      
      Comp_1 <- sapply(CppSim$simulationResults, function(x){x[[cmpName]][,1]})
      Comp_2 <- sapply(CppSim$simulationResults, function(x){x[[cmpName]][,2]})
      Comp_3 <- sapply(CppSim$simulationResults, function(x){x[[cmpName]][,3]})
      Comp_4 <- sapply(CppSim$simulationResults, function(x){x[[cmpName]][,4]})
      Comp1mean <- apply(Comp_1, 1, median)
      Comp2mean <- apply(Comp_2, 1, median)
      Comp3mean <- apply(Comp_3, 1, median)
      Comp4mean <- apply(Comp_4, 1, median)
      Comp1mean_R <- apply(Comp_1_R, 1, median)
      Comp2mean_R <- apply(Comp_2_R, 1, median)
      Comp3mean_R <- apply(Comp_3_R, 1, median)
      Comp4mean_R <- apply(Comp_4_R, 1, median)
      
      Comp1LB <- apply(Comp_1, 1, quantile, probs = c(0.1))
      Comp2LB <- apply(Comp_2, 1, quantile, probs = c(0.1))
      Comp3LB <- apply(Comp_3, 1, quantile, probs = c(0.1))
      Comp4LB <- apply(Comp_4, 1, quantile, probs = c(0.1))
      
      Comp1UB <- apply(Comp_1, 1, quantile, probs = c(0.9))
      Comp2UB <- apply(Comp_2, 1, quantile, probs = c(0.9))
      Comp3UB <- apply(Comp_3, 1, quantile, probs = c(0.9))
      Comp4UB <- apply(Comp_4, 1, quantile, probs = c(0.9))
      
      
      Comp1LB_R <- apply(Comp_1_R, 1, quantile, probs = c(0.1))
      Comp2LB_R <- apply(Comp_2_R, 1, quantile, probs = c(0.1))
      Comp3LB_R <- apply(Comp_3_R, 1, quantile, probs = c(0.1))
      Comp4LB_R <- apply(Comp_4_R, 1, quantile, probs = c(0.1))
      
      Comp1UB_R <- apply(Comp_1_R, 1, quantile, probs = c(0.9))
      Comp2UB_R <- apply(Comp_2_R, 1, quantile, probs = c(0.9))
      Comp3UB_R <- apply(Comp_3_R, 1, quantile, probs = c(0.9))
      Comp4UB_R <- apply(Comp_4_R, 1, quantile, probs = c(0.9))
      
      ckfc <- function(a,b){
        passed <- (var(a) == var(b) ||
        cor(a, b) >= 0.95 &&
        mean(abs(a-b))/mean(b) < 0.1)
        expect(passed, "R and C++ simulations do not agree.")
        passed
      }
      
      if (!(ckfc(Comp1LB, Comp1LB_R) &&
            ckfc(Comp2LB, Comp2LB_R) &&
            ckfc(Comp3LB, Comp3LB_R) &&
            ckfc(Comp4LB, Comp4LB_R) &&
            ckfc(Comp1UB, Comp1UB_R) &&
            ckfc(Comp2UB, Comp2UB_R) &&
            ckfc(Comp3UB, Comp3UB_R) &&
            ckfc(Comp4UB, Comp4UB_R) &&
            ckfc(Comp1mean, Comp1mean_R) &&
            ckfc(Comp2mean, Comp2mean_R) &&
            ckfc(Comp3mean, Comp3mean_R) &&
            ckfc(Comp4mean, Comp4mean_R))){
        warning("Simulations don't match!")
      }
    }
    compareComp(cmpName = "I", Rcomp, CppSim)
})

