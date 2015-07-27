# dataModel module helper function
DataModel = function(Y, type = c("identity", "overdispersion"), compartment = c("I_star", "R_star"), phi = NA)
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
    structure(list("Y"=Y, 
                   "type"=type,
                   "compartment"=compartment,
                   "phi"=phi), class = "DataModel")
}

