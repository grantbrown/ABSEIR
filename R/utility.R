Ifelse <- function(x, a, b)
{
    if (x){
        return(a)
    }
    b
}

validateIf <- function(cond, validator)
{
    if (cond){
        return(validator)
    }
    return(function(x){})
}

mustHaveMember <- function(name){
    function(arg){
        if (!(name %in% names(arg)))
        {
            return(paste("must have element: ", name, sep = ""))
        }
    }
}


mustBeOneOf <- function(elements){
    function(arg){
        if (!(arg[1] %in% elements))
        {
            return(paste("must be one of: ", paste0(elements, collapse = ", "), sep = ""))
        }
    }
}

mustBeLen <- function(len)
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

mustHaveClass <- function(classes)
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

mustBeInRange <- function(lower=-Inf,upper=Inf, 
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
                          ...)
{
    argValidators <- list(...)
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
}
