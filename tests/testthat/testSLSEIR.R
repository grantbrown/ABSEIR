test_that("Single location SEIR models can be constructed", {
  data(Kikwit1995)
  data_model_1 = DataModel(Kikwit1995$Count,
                         type = "identity",
                         compartment="I_star",
                         cumulative=FALSE)
  data_model_2 = DataModel(cumsum(Kikwit1995$Count),
                         type = "identity",
                         compartment="I_star",
                         cumulative=TRUE)
  data_model_3 = DataModel(Kikwit1995$Count,
                           type = "overdispersion",
                           compartment="I_star",
                           cumulative=FALSE, params=list(phi = 1))
  dataModelList = list(data_model_1,
                       data_model_2,
                       data_model_3)
  
  intervention_term = cumsum(Kikwit1995$Date >  as.Date("05-09-1995", "%m-%d-%Y"))
  intervention_term = intervention_term/max(intervention_term)
  exposure_model = ExposureModel(cbind(1,intervention_term),
                                   nTpt = nrow(Kikwit1995),
                                   nLoc = 1,
                                   betaPriorPrecision = 0.5,
                                   betaPriorMean = 0)
  # There's no reinfection in this case, so we just use a "SEIR" model. 
  reinfection_model = ReinfectionModel("SEIR")
  # we have no distance model, because it's a non-spatial analysis, but we still
  # need to create an empty DistanceModel for ABSEIR. This approach is a bit clunky, 
  # so we're open to suggestions. One option may be the implmentation of new verisons
  # of the qSEIR functions included in libSpatialSEIR. 
  distance_model = DistanceModel(list(matrix(0)))
  
  # Set initial population sizes
  initial_value_container = InitialValueContainer(S0=5.36e6,
                                                  E0=2,
                                                  I0=2,
                                                  R0=0)
  
  
  # Model to describe E to I and I to R transition probabilities. Prior values from 
  # Lekone and Finkenstadt 2006
  # This is going to be expanded soon, so that we won't need to rely on 
  # the exponential assumption
  transition_priors_1 = ExponentialTransitionPriors(p_ei = 1-exp(-1/5), 
                                                  p_ir= 1-exp(-1/7),
                                                  p_ei_ess = 100,
                                                  p_ir_ess = 100)
  transition_priors_2 = PathSpecificTransitionPriors(Z1 = function(x){dunif(x,3,7)},
                                                     Z2 = function(x){dunif(x, 7, 28)})
  transition_priors_3 = WeibullTransitionPriors(10,10,10,10,10,10,10,10)
  
  transitionPriorsList = list(transition_priors_1, transition_priors_2, transition_priors_3)
  # Set algorithm configuration
  sampling_control = SamplingControl(seed = 123123, 
                                     n_cores = 2,
                                     algorithm="Beaumont2009",
                                     list(batch_size = 100,
                                          epochs = 5,
                                          max_batches = 2,
                                          shrinkage = 0.99,
                                          multivariate_perturbation=FALSE
                                     )
  )
  
  results <- as.list(1:(length(dataModelList)*length(transitionPriorsList)))
  i <- 1
  for (dm in 1:length(dataModelList)){
    for (tp in 1:length(transitionPriorsList))
    {
      tryCatch({
        result = SpatialSEIRModel(dataModelList[[dm]],
                                   exposure_model,
                                   reinfection_model,
                                   distance_model,
                                   transitionPriorsList[[tp]],
                                   initial_value_container,
                                   sampling_control,
                                   samples = 100,
                                   verbose = FALSE)
      
        simulated = epidemic.simulations(result, replicates = 25)
        i <- i + 1
        results[[i]] = list(results=result, simulationss <- simulated)
      }, warning = function(w){
        print(w)
        expect_equal(paste(
          "Warnings generated while fitting single location model: (", 
          dm, ", ", "tp", ")\n", sep = ""), FALSE)}, 
      error = function(e){
        print(e)
        expect_equal(paste("Errors generated while fitting single location model: (", 
                     dm, ", ", "tp", ")\n", sep = ""), FALSE)
      }
      )
    }
  }
})
