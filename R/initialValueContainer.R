#' Create an InitialValueContainer object, which catalogues the starting population
#' values for an epidemic.
#'
#' @param S0 An integer vector of susceptible counts
#' @param E0 An integer vector of exposed counts
#' @param I0 An integer vector of infectious counts
#' @param R0 AN integer vector of removed counts
#' @details
#'  While approximate Bayesian techniques avoid the need to specify starting values
#'  for the spatiotemporally varying population values (i.e. the current number of 
#'  susceptible, exposed, infectious, and removed individuals), we still need
#'  to specify how many individuals start out in each of these categories. In the 
#'  future, we may support the specification of priors on these terms, but for
#'  now they are considered fixed. 
#'
#' @examples inits <- InitialValueContainer(c(1000,1000),
#'                                           c(10,0),
#'                                           c(0,0),
#'                                           c(100,100))
#' @export
InitialValueContainer = function(S0, E0, I0, R0) 
{
    # get rid of any data frame nonsense
    structure(list(S0 = as.numeric(S0),
                   E0 = as.numeric(E0),
                   I0 = as.numeric(I0),
                   R0 = as.numeric(R0)), class="InitialValueContainer")
}



