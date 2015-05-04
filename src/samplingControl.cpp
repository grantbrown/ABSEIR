#include <Rcpp.h>
#include <samplingControl.hpp>


using namespace Rcpp;

samplingControl::samplingControl(SEXP width, SEXP seed, SEXP cores)
{
    Rcpp::IntegerVector in_width(width); 
    Rcpp::IntegerVector in_seed(seed);
    Rcpp::IntegerVector in_cores(cores);
    if (in_width.length() != 1)
    {
        ::Rf_error("Simulation width must be of length 1.");
    }
    if (in_seed.length() != 1)
    {
        ::Rf_error("Simulation seed must be of length 1.");
    }
    if (in_cores.length() != 1)
    {
        ::Rf_error("Number of cores must be of length 1.");
    }

    simulation_width = in_width[0];
    random_seed = in_seed[0];
    CPU_cores = in_cores[0];
}

samplingControl::~samplingControl()
{
    if (prot !=0 ){
        ::Rf_error("can't delete samplingControl, still being used.\n");
    }
}

int samplingControl::getModelComponentType()
{
    return(LSS_SAMPLING_CONTROL_MODEL_TYPE);
}

RCPP_MODULE(mod_samplingControl)
{
    using namespace Rcpp;
    class_<samplingControl>( "samplingControl" )
    .constructor<SEXP, SEXP, SEXP>();
}


