#' Plot a graphical summary of the marginal posterior distribution of SEIR model parameters
#' @param x a \code{\link{SpatialSEIRModel}} object
#' @param \dots additional arguments to be passed to plotting functions. 
#' @examples \dontrun{plot(modelObject)}
#' @export
#' @import stats
#' @importFrom graphics  "hist"
plot.SpatialSEIRModel = function(x, ...)
{
    for (i in 1:ncol(x$param.samples))
    {
        cat ("Press return to see next plot:")
        hist(x$param.samples[,i], 
             main = colnames(x$param.samples)[i],
             xlab = "Parameter Value",
             ylab = "MPD",
             freq=FALSE,
             ...
             )
        if (i != ncol(x$param.samples)) line <- readline()
    }
}

