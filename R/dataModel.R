# dataModel module helper function
DataModel = function(Y, type = c("identity", "overdispersion"), compartment = c("I_star", "R_star"), phi = NA)
{
    if (compartment != "I_star")
    {
        stop(paste("Currently, only the I_star compartment is supported.", 
             " Alternate models are planned, so if you find yourself",
             " in need of such features, please email grant-brown@uiowa.edu"))
    }
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
    na_mask = is.na(Y)
    Y[na_mask] = -Inf
    structure(list("Y"=Y, 
                   "type"=type,
                   "compartment"=compartment,
                   "phi"=phi,
                   "na_mask" = na_mask), class = "DataModel")
}

