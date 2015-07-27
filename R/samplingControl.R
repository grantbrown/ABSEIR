# samplingControl module helper function
SamplingControl = function(sim_width, seed, n_cores)
{
    structure(list("sim_width" = sim_width,
                   "seed" = seed,
                   "n_cores" = n_cores), class = "SamplingControl")
}


