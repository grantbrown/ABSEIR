Ifelse <- function(x, a, b)
{
    if (x){
        return(a)
    }
    b
}

argLengthValidator <- function(len)
{
    function(arg)
    {
        err <- list()
        if (length(arg) < min(len) || length(arg) > max(len))
        {
            return(paste("invalid argument length: ", length(arg), sep = ""))
        }
    }
}

argClassValidator <- function(classes)
{
    function(arg){
        if (length(intersect(class(arg), classes)) == 0) 
        {
            return(paste("Object not of class: ", 
                         paste(classes, collapse = ", "),
                         ", is: ", paste(class(arg), collapse = ", "), 
                         sep = ""))
        }
    }
}

numericRangeValidator <- function(lower=-Inf,upper=Inf, 
                                  rightClosed = TRUE, 
                                  leftClosed=TRUE){
    rightComp <- ifelse(rightClosed, `<=`, `<`)
    leftComp <- ifelse(leftClosed, `>=`, `>`)

    function(arg){
        err = list()
        tryCatch({
            if ((!leftComp(arg, lower)) || (!rightComp(arg, upper)))
            {
                err[[length(err)+1]] <- paste("Range Error: must be in (", lower, 
                                              ", ", upper, ")", sep = "")
            }  
            return(err)
            }, 
            warning = function(w){
                err[[length(err)+1]] <- w
            },
            error = function(e){
                err[[length(err)+1]] <- e
            }
        )
        return(err)
    }
}

checkArgument <- function(argName, 
                          argValidators,
                          dotValidators,
                          ...)
{
    argValue <- get(argName, parent.frame())
    errs <- 0;
    for (validator in argValidators)
    {
       validationResult <- validator(argValue) 
       if (length(validationResult) > 0)
       {
           errs <- errs + 1
           cat(paste("Argument validation error for:", argName, ".\n"))
           for (r in validationResult){
               print(validationResult)
           }
       }
    }
    if (errs > 0)
    {
        stop("Invalid argument.")
    }
    if (!missing(dotValidators))
    {
        warning("checkArgument does not currently support dotValidators")
    }
}
