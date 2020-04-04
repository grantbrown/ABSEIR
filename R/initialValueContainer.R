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
InitialValueContainer = function(S0, E0, I0, R0, 
                                 type = c("identity", "uniform"),
                                 params = list()) 
{
    type <- type[1]
    checkArgument("type", mustBeOneOf(c("identity", "uniform")))
    checkArgument("params", 
                  mustHaveClass("list"),
                  validateIf(type == "uniform", 
                             mustHaveMember("max_S0")),
                  validateIf(type == "uniform", 
                             mustHaveMember("max_E0")),
                  validateIf(type == "uniform", 
                             mustHaveMember("max_I0")),
                  validateIf(type == "uniform", 
                             mustHaveMember("max_R0"))
    )
    l1 <- length(as.numeric(S0))
    checkArgument("E0", mustBeLen(l1))
    checkArgument("I0", mustBeLen(l1))
    checkArgument("R0", mustBeLen(l1))
    
  
    if (type == "identity"){
    
      
      return(# get rid of any data frame nonsense
        structure(list(S0 = as.numeric(S0),
                     E0 = as.numeric(E0),
                     I0 = as.numeric(I0),
                     R0 = as.numeric(R0), 
                     max_S0 = as.numeric(S0), 
                     max_E0 = as.numeric(E0),
                     max_I0 = as.numeric(I0),
                     max_R0 = as.numeric(R0),
                     type = 1), class="InitialValueContainer")
      
      )
    } else if (type == "uniform"){
      return(# get rid of any data frame nonsense
        structure(list(S0 = as.numeric(S0),
                       E0 = as.numeric(E0),
                       I0 = as.numeric(I0),
                       R0 = as.numeric(R0), 
                       max_S0 = params$max_S0, 
                       max_E0 = params$max_E0,
                       max_I0 = params$max_I0,
                       max_R0 = params$max_R0, 
                       type = 2), class="InitialValueContainer")
        
      )
    } else{
      stop("Bad user, bad!") # shouldn't actually be possible to hit this. 
  }
}



