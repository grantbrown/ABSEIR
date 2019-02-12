test_that("Epsilon calculations are correct",{
    library(ABSEIR)
    data(Kikwit1995)
  
    data_model_1 = DataModel(Kikwit1995$Count,
                           type = "identity",
                           compartment="I_star",
                           cumulative=FALSE)
    
    data_model_2 = DataModel(Kikwit1995$Count,
                             type = "overdispersion",
                             compartment="I_star",
                             cumulative=FALSE,
                             params = list(phi = 1))
    
    data_model_3 = DataModel(cumsum(Kikwit1995$Count),
                             type = "identity",
                             compartment="I_star",
                             cumulative=TRUE)
    
    data_model_4 = DataModel(cumsum(Kikwit1995$Count),
                             type = "overdispersion",
                             compartment="I_star",
                             cumulative=TRUE,
                             params = list(phi = 1))
    data_model_5 = DataModel(Kikwit1995$Count,
                             type = "fractional",
                             compartment="I_star",
                             cumulative=FALSE,
                             params = list(report_fraction = 0.9,
                                           report_fraction_ess=100))
    data_model_6 = DataModel(cumsum(Kikwit1995$Count),
                             type = "fractional",
                             compartment="I_star",
                             cumulative=TRUE,
                             params = list(report_fraction = 0.9,
                                           report_fraction_ess=100))
    
    
    dataModelList <- list(data_model_1,
                          data_model_2,
                          data_model_3,
                          data_model_4,
                          data_model_5,
                          data_model_6)
    
    intervention_term = cumsum(Kikwit1995$Date >  as.Date("05-09-1995", "%m-%d-%Y"))
    exposure_model = ExposureModel(cbind(1,intervention_term),
                                     nTpt = nrow(Kikwit1995),
                                     nLoc = 1,
                                     betaPriorPrecision = 0.5,
                                     betaPriorMean = 0)
    
    
    # There's no reinfection in this case, so we just use a "SEIR" model. 
    reinfection_model = ReinfectionModel("SEIR")
    
    distance_model = DistanceModel(list(matrix(0)))
    
    initial_value_container = InitialValueContainer(S0=5.36e6,
                                                    E0=2,
                                                    I0=2,
                                                    R0=0)
    
    transition_priors = ExponentialTransitionPriors(p_ei = 1-exp(-1/5), 
                                                    p_ir= 1-exp(-1/7),
                                                    p_ei_ess = 10,
                                                    p_ir_ess = 10)
    
    # Set algorithm configuration
    sampling_control = SamplingControl(seed = 123123, 
                                       n_cores = 8,
                                       algorithm="Beaumont2009",
                                       list(batch_size = 2000,
                                            epochs = 2,
                                            max_batches = 100,
                                            shrinkage = 0.99,
                                            multivariate_perturbation=FALSE,
                                            keep_compartments = TRUE
                                       )
    )
    
    rf <- function(dmobj){
      runtime = system.time(result <- SpatialSEIRModel(dmobj,
                                                         exposure_model,
                                                         reinfection_model,
                                                         distance_model,
                                                         transition_priors,
                                                         initial_value_container,
                                                         sampling_control,
                                                         samples = 100,
                                                         verbose = 2))
     e1 <- result$epsilon
     e2 <- sapply(result$simulationResults, function(x){x$result})
     expect(all(e1 == e2), "Epsilon calculations don't agree between simulation and inference code")
     if (dmobj$cumulative != TRUE){
       e3 <- sapply(result$simulationResults, function(x){
         sum((x$I_star - Kikwit1995$Count)^2)
       })
     }
     else{
       e3 <- sapply(result$simulationResults, function(x){
         sum((cumsum(x$I_star) - cumsum(Kikwit1995$Count))^2)
       })
     }
     
     # 554 621 595 554 574 638 806 575 565 467 
     if (dmobj$type == "identity"){
      expect(all(e1 == e3), "Epsilon calculations don't match squared-distance")
     } else{
       expect(!all(e1 == e3), "Epsilon calculations unexpectedly match squared-distance")
     }
    }
    lapply(dataModelList, rf)
})
