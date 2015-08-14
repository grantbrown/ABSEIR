\name{SpatialSEIRModel}
\alias{SpatialSEIRModel}
\title{Fit a single locaton or spatial SEIR model using ABC rejection sampling}
\description{This is the main analysis function provided by \code{ABSEIR}, 
        combining each of the model components created by 
        \code{DataModel}, \code{ExposureModel}, \code{ReinfectionModel},
        \code{TransitionPriors}, \code{InitialValueContainer}, and 
        \code{SamplingControl}, and accepting additional fitting arguments
        to produce samples from the approximate posterior distribution.}
\usage{
    SpatialSEIRModel(data_model, exposure_model, reinfection_model,
        distance_model, transition_priors, initial_value_container,
        sampling_control, samples = 100, verbose = FALSE)
}
\arguments{
  \item{data_model}{A data model object, describing the link between
  the observed data and the unobserved epidemic counts. Valid data models
  are created using the \code{\link{DataModel}} function.}
  \item{exposure_model}{An exposure model object, which describes the
  spatial and temporal variability of the exposure/infection process. Valid 
  exposure models are created using the \code{\link{ExposureModel}} function.}
  \item{reinfection_model}{A reinfection model object, which describes
  whether or not individuals are able to return from the Removed category
  to the Susceptible population. Valid reinfection models are created using the
  \code{\link{ReinfectionModel}} function.}
  \item{distance_model}{A distance model object which describes the underlying
  contact network, in addition to prior parameters which constrain the contact
  process. Valid distance models are created using \code{\link{DistanceModel}}}
  \item{transition_priors}{An object containing information about the E to I and
  I to R transition prior parameters. These are created using the 
  \code{\link{TransitionPriors}} function by specifying transition
  probabilities and corresponding effective sample sizes (amount of prior 
  information).} 
  \item{initial_value_container}{An object specifying the initial state
  of the epidemic for each spatial location, created by the
  \code{\link{InitialValueContainer}} function.} 
  \item{sampling_control}{An object specifying information about the sampling
  algorithm. In particular, the sampling_control argument should specify the
  number of CPU cores to employ, and the random seed to use. 
  Sampling control objects are created by the 
  \code{\link{SamplingControl}} function.}
  \item{samples}{This argument determines the number of posterior samples
  to be returned.}
  \item{verbose}{Logical - \code{verbose==TRUE} indicates that ABSEIR should
  provide diagnostic and convergence information.}
}
\details{
    Numerous ABC algorithms have been developed, but as of now ABSEIR provides 
    just two. The first algorithm is the basic rejection algorithm of Rubin 
    1980. While this approach performs well when good prior information is available, 
    it can be extremely inefficient when prior distributions are diffuse with 
    respect to the posterior. To address this shortcoming, we have implemented
    the Sequential Monte-Carlo approach proposed by Beaumont 2009, 2010. We may
    provide additional algorithms in the future, in particular that of Del Moral
    et al. 2012. 
}
\value{
  An object with S3 class \code{"SpatialSEIRModel"} containing:
  \item{param.samples}{A matrix containing accepted samples from the prior 
  distribution.}
  \item{epsilon}{Relative similarity value between simulation and 
  observed data. 
  }
  \item{modelComponents}{A copy of the model component arguments accepted as 
  arguments to \code{SpatialSEIRModel} to enable additional sampling or 
  compartment simulation.}
}
\references{
    Beaumont, M.A., Cornuet, J.-M., Marin, J.-M. et al., 2009. 
    Adaptive approximate Bayesian computation. Biometrika 96,983-990.

    Del Moral, P., Doucet, A., Jasra, A., 2012. An adaptive sequential Monte Carlo
    method for approximate Bayesian computation. Statistics and Computing 22, 1009-1020.

    Rubin, D., 1980. Bayesianly justifiable and relevant frequency calculations 
    for the applied statistician. The annals of Statistics 12, 1151-1172.         
}

\author{Grant Brown <grant-brown@uiowa.edu>}
\seealso{\code{\link{DataModel}}, \code{\link{ExposureModel}},
         \code{\link{ReinfectionModel}}, 
         \code{\link{InitialValueContainer}},
         \code{\link{DistanceModel}}, 
         \code{\link{TransitionPriors}}
         }
\examples{
    print("Somebody tell Grant to make some examples.")
}
