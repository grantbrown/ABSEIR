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
        sampling_control, samples = 100, accept_fraction = 1e-4,
        batch_size = 50000, verbose = FALSE)
}
\arguments{
  \item{data_model}{A data model object, describing the link between
  the observed data and the unobserved epidemic counts. Valid data models
  are created using the \code{\link{DataModel}} function.}
  \item{exposure_model}{An exposure model object, which describes the
  spatial and temporal variability of the exposure/infection process. Valid 
  exposure models are created using the \code{\link{ExposureModel}}} function.  
  \item{reinfection_model}{A reinfection model object, which describes
  whether or not individuals are able to return from the Removed category
  to the Susceptible population. Valid reinfection models are created using the
  \code{\link} ReinfectionModel function.}
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
  number of CPU cores to employ. Sampling control objects are created by the 
  \code{\link{SamplingControl}} function.}
  \item{samples}{This argument determines the number of posterior samples
  to be returned.}
  \item{accept_fraction}{Numeric - The accept_fraction argument determines what 
  proportion of proposed parameters are accepted. An acceptance fraction of 
  1 produces samples from the prior distribution. Good values of this rate
  are problem specific; starting with a relatively liberal value like 1e-2 
  and moving towards zero as computation time permits generally works well.
  ABSEIR will simulate \code{ceil(samples/accept_fraction)} total epidemics. 
  }
  \item{batch_size}{Integer - For systems in which memory is an issue, ABSEIR does not
  generate all required parameter combinations at once. Instead, the package
  generates \code{batch_size} samples from the prior distribution at a time, and 
  keeps only \code{samples} of them. Subsequent required batches are merged into
  this initial sample.}
  \item{verbose}{Logical - \code{verbose==TRUE} indicates that ABSEIR should
  provide diagnostic information.}
}
\details{
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
\references{Manuscript in progress}

\author{Grant Brown <grant-brown@uiowa.edu>}
\seealso{\code{\link{DataModel}}, \code{\link{ExposureModel}}}
\examples{
    print("examples go here")
}
