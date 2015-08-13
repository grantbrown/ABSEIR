# samplingControl module helper function
SamplingControl = function(seed, n_cores, algorithm="Beaumont2009",
                           params=NA)                           
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

    # Set default parameters where needed
    if (class(params) == "logical" && all(is.na(params))){
        if (algorithm == "Beaumont2009")
        { 
            params = list(acceptance_fraction = -1,
                 batch_size = 5000,
                 epochs = 100, 
                 shrinkage = 0.9,
                 max_batches = 20)
        }
        else
        {
            params = list(acceptance_fraction = 0.01,
                 batch_size = 10000,
                 epochs = 1,
                 shrinkage = 1,
                 max_batches = 1)
        }
    }
    else if (class(params) == "list")
    {
        if (algorithm == "Beaumont2009")
        {
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 5000
            }
            if (!("epochs" %in% names(params))) {
                params[["epochs"]] = 100
            }
            if (!("shrinkage" %in% names(params))) {
                params[["shrinkage"]] = 0.9
            }
            if (!("acceptance_fraction" %in% names(params))) {
                params[["acceptance_fraction"]] = 1
            }
            if (!("max_batches" %in% names(params))) {
                params[["max_batches"]] = 20
            }
        }
        else if (algorithm == "BasicABC")
        {
            if (!("batch_size" %in% names(params))) {
                params[["batch_size"]] = 10000
            }
            if (!("epochs" %in% names(params))) {
                params[["epochs"]] = 1
            }
            if (!("shrinkage" %in% names(params))) {
                params[["shrinkage"]] = 1
            }
            if (!("acceptance_fraction" %in% names(params))) {
                params[["acceptance_fraction"]] = 0.01
            }
            if (!("max_batches" %in% names(params))) {
                params[["max_batches"]] = 1
            }
        }
        else
        {
            stop("params must be a list.")
        }
    }


    structure(list("sim_width" = 1,
                   "seed" = seed,
                   "n_cores" = n_cores,
                   "acceptance_fraction" = params$acceptance_fraction,
                   "batch_size" = params$batch_size,
                   "algorithm" = alg,
                   "epochs" = params$epochs,
                   "shrinkage" = params$shrinkage,
                   "max_batches" = params$max_batches
                   ), class = "SamplingControl")
}


