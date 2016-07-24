#include <Rcpp.h>
#include <Eigen/Core>
#include <RcppEigen.h>
#include <cmath>
#include <math.h>
#include <spatialSEIRModel.hpp>
#include <dataModel.hpp>
#include <exposureModel.hpp>
#include <reinfectionModel.hpp>
#include <distanceModel.hpp>
#include <transitionPriors.hpp>
#include <initialValueContainer.hpp>
#include <samplingControl.hpp>
#include <util.hpp>
#include <SEIRSimNodes.hpp>

Rcpp::List spatialSEIRModel::sample_basic(int nSample, int verbose, 
                                          std::string sim_type_atom)
{
    Rcpp::stop("Basic rejection algorithm not currently supported.\n");
    // Return dummy value
    Rcpp::List outList;
    return(outList);
}
