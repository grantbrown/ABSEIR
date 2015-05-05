## distanceModel module helper function
buildDistanceModel = function(distanceList, 
                              scaleMode = c("none","rowscale","invsqrt"),
                              priorAlpha=1.0,
                              priorBeta=1.0)
{
    scaleMode = scaleMode[1]
    rowScale = function(mat)
    {
        mat/matrix(apply(mat,1,sum), nrow = nrow(mat), ncol = ncol(mat))
    }

    invSqrt = function(mat)
    {
        matrix(ifelse(mat == 0, 0, 1/sqrt(mat)), nrow = nrow(mat), ncol = ncol(mat))
    }


    # Check if we've got a single matrix
    if (class(distanceList) == "matrix")
    {
        distanceList = list(distanceList)
    }
    # Make sure that at this point we have a list of matrices. 

    if (class(distanceList) != "list")
    {
        stop("Error: distanceList must be a list of matrices.")
    }

    # Check for valid matrices
    distanceDim = NA
    for (i in 1:length(distanceList))
    {
        if (class(distanceList[[i]]) != "matrix")
        {
            stop("Distance metrics must be matrices.")
        }
        distanceDim = dim(distanceList[[1]])
        newDim = dim(distanceList[[i]])
        if (any(newDim != distanceDim) || (newDim[1] != newDim[2]))
        {
            stop("Distance matrices must be square and of the same dimension.")
        }
        if (any(is.na(distanceList[[i]])))
        {
            stop("NA's are not allowed in distance matrices.")
        }
        if (scaleMode == "rowscale")
        {
            distanceList[[i]] = rowScale(distanceList[[i]])
        }
        if (scaleMode == "invsqrt")
        {
            distanceList[[i]] = invSqrt(distanceList[[i]])
        }
    } 

    finalDistanceModel = new( distanceModel )
    for (i in 1:length(distanceList))
    {
        finalDistanceModel$addDistanceMatrix(distanceList[[i]])
    }
    finalDistanceModel$setPriorParameters(priorAlpha, priorBeta);
    finalDistanceModel
}

# dataModel module helper function
buildDataModel = function(Y, type = c("identity", "overdispersion"), compartment = c("I_star", "R_star"), phi = NA)
{
    if (class(Y) != "matrix")
    {
        Y = as.matrix(Y)
    }
    type = type[1]
    if (length(phi) == 1 && is.na(phi) && type != "identity")
    {
        stop("Must specify overdispersion parameter (phi) for non-identity data model.")
    }
    else if (length(phi) == 1 && is.na(phi))
    {
        return(new(dataModel, Y, type, compartment, -1.0))
    }
    else if (class(phi) != "numeric" || length(phi) != 1)
    {
        stop("Non identy data model currently requires a single overdispersion parameter") 
    }
    outModel = new(dataModel, Y, type, compartment, phi)
    outModel
}

# reinfectionModel module helper function
buildReinfectionModel = function(reinfectMode = c("SEIR", "SEIRS", "Fixed"),X_prs = NA, 
                                 betaPrs=NA , priorPrecision = NA, priorMean = NA)
{
    integerMode = ifelse(reinfectMode[1] == "SEIR", 3, 
                  ifelse(reinfectMode[1] == "SEIRS", 1, 
                  ifelse(reinfectMode[1] == "Fixed", 2, NA)))
    if (is.na(integerMode))
    {
        stop(paste("Invalid mode: ", reinfectMode[1], sep = ""))
    }
    if (integerMode != 3 && (is.na(X_prs) || is.na(betaPrs)))
    {
        stop("If reinfection mode is not SEIR, X_prs and betaPrs must be supplied.")
    }
    if (integerMode == 1 && (is.na(priorPrecision) || is.na(priorMean)))
    {
        stop("If reinfection parameters are going to be estimated, priorPrecision and priorMean must be specified.")
    }

    reinfectionmod = new(reinfectionModel, integerMode);
    if (integerMode != 3)
    {
        if (class(X_prs) != "matrix")
        {
            X_prs = as.matrix(X_prs)
        }
        if (all(is.na(priorPrecision)))
        {
           priorPrecision = rep(0.1, ncol(X_prs))     
        }
        else if (length(priorPrecision) == 1)
        {
            priorPrecision = rep(priorPrecision, ncol(X_prs))               
        }
        if (all(is.na(priorMean)))
        {
            priorMean = rep(0, ncol(X_prs))
        }
        else if (length(priorMean) == 1)
        {
            priorMean = rep(priorMean, ncol(X_prs))
        }
        reinfectionmod$buildReinfectionModel(X_prs, betaPrs, priorMean, priorPrecision);
    }
    reinfectionmod
}

# samplingControl module helper function
buildSamplingControl = function(sim_width, seed, n_cores)
{
    samplingControlInstance = new ( samplingControl, sim_width, seed, n_cores)
    samplingControlInstance        
}

# transitionPriors module helper functions
buildTransitionPriorsFromProbabilities = function(p_ei, p_ir, p_ei_ess, p_ir_ess)
{
    tp = new(transitionPriors)
    tp$setPriorsFromProbabilities(p_ei, p_ir, p_ei_ess, p_ir_ess)
    tp
}
buildTransitionPriorsManually = function(priorAlpha_gammaEI, priorBeta_gammaEI,
                                         priorAlpha_gammaIR, priorBeta_gammaIR)
{
     tp = new(transitionPriors)
     tp$setPriorsManually(priorAlpha_gammaEI, priorBeta_gammaEI,
                          priorAlpha_gammaIR, priorBeta_gammaIR)
     tp
}
buildUniformTransitionPriors = function()
{
    tp = new(transitionPriors)
    tp$setUniformPriors()
    tp
}

# exposureModel module helper function
buildExposureModel = function(X,nTpt, nLoc,betaPriorPrecision=NA,
                              betaPriorMean=NA,offset=NA)
{
    nBeta = ncol(X)
    if (class(X) != "matrix")
    {
        print("Warning: X should be a matrix.")
    }
    if (nrow(X) != nTpt * nLoc){
        stop("Invalid data dimensions: X should be a matrix composed of nLoc row-wise blocks of dimension nTpt*p.") 
    } 
    if (length(betaPriorPrecision) == 1 && is.na(betaPriorPrecision))
    {
        print("No prior precision specified, using 0.1")
        betaPriorPrecision = rep(0.1, nBeta)
    }
    else if (length(betaPriorPrecision) == 1)
    {
        betaPriorPrecision = rep(betaPriorPrecision, nBeta)
    }
    if (length(betaPriorMean) == 1 && is.na(betaPriorMean))
    {
        print("No prior mean specified, using zero.")
        betaPriorMean = rep(0, nBeta)
    }
    else if (length(betaPriorMean) == 1)
    {
        betaPriorMean = rep(betaPriorMean, nBeta)
    }
    if (length(offset) == 1 && is.na(offset))
    {
        print("Assuming equally spaced count data.")
    }
    else if (length(offset) == 1){
        offset = rep(offset, nTpt)
    }
    else if (length(offset) != nTpt){
        stop("Offset must be the of length nTpt")
    }
    ExposureModel = new(exposureModel,X,nTpt,nLoc,betaPriorMean,betaPriorPrecision)
    if (all(!is.na(offset)))
    {
        ExposureModel$offsets = offset
    }
    ExposureModel
}

# initialValueContainer module helper function
buildInitialValueContainer = function(S0, E0, I0, R0) 
{
    # get rid of any data frame nonsense
    S0 = as.numeric(S0)
    E0 = as.numeric(E0)
    I0 = as.numeric(I0)
    R0 = as.numeric(R0)

    InitialValueContainer = new(initialValueContainer)
    InitialValueContainer$setInitialValues(S0, E0, I0, R0)
    InitialValueContainer
}


# SEIRModel module helper function
buildSEIRModel = function(dataModelInstance,
                          exposureModelInstance,
                          reinfectionModelInstance,
                          distanceModelInstance,
                          transitionPriorsInstance,
                          initialValueContainer,
                          samplingControlInstance)
{

    interface = new( spatialSEIRModel, dataModelInstance,
                                   exposureModelInstance,
                                   reinfectionModelInstance,
                                   distanceModelInstance,
                                   transitionPriorsInstance,
                                   initialValueContainer,
                                   samplingControlInstance)
    interface
}
