# samplingControl module helper function
SamplingControl = function(sim_width, seed, n_cores, 
                           acceptance_fraction,
                           batch_size, algorithm="BasicABC",
                           epochs = 100)
{
    alg = ifelse(algorithm == "BasicABC", 1,
          ifelse(algorithm == "Beaumont2009", 2, 
                 NA))
    if (is.na(alg))
    {
        stop(paste("Only the basic rejection algorithm and a modified version",
                   "of the SMC algorithm from Beaumont 2009 are currently",  
                   "supported.\n", sep = ""))
    }
    structure(list("sim_width" = sim_width,
                   "seed" = seed,
                   "n_cores" = n_cores,
                   "acceptance_fraction" = acceptance_fraction,
                   "batch_size" = batch_size,
                   "algorithm" = alg,
                   "epochs" = epochs), class = "SamplingControl")
}


