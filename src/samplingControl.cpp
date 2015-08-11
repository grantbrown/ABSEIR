#include <Rcpp.h>
#include <samplingControl.hpp>


using namespace Rcpp;

samplingControl::samplingControl(SEXP width, SEXP seed, SEXP cores, 
                                 SEXP frac, SEXP batch, SEXP alg,
                                 SEXP n_epochs)
{
    Rcpp::IntegerVector in_width(width); 
    Rcpp::IntegerVector in_seed(seed);
    Rcpp::IntegerVector in_cores(cores);
    Rcpp::IntegerVector in_algorithm(alg);
    Rcpp::NumericVector in_frac(frac);
    Rcpp::IntegerVector in_batch(batch);
    Rcpp::IntegerVector in_epochs(n_epochs);

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
    if (in_algorithm.length() != 1 || (in_algorithm(0) != ALG_BasicABC && 
                in_algorithm(0) != ALG_ModifiedBeaumont2009))
    {
        ::Rf_error("Algorithm specification must be of length 1 and equal to 1 or 2.");
    }
    if (in_batch.length() != 1)
    {
        ::Rf_error("Batch size must be of length 1.");
    }
    if (in_frac.length() != 1)
    {
        ::Rf_error("Acceptance fraction must be of length 1.");
    }
    if (in_epochs.length() != 1)
    {
        ::Rf_error("Epochs must be of length 1.");
    }

    simulation_width = in_width(0);
    random_seed = in_seed(0);
    CPU_cores = in_cores(0);
    algorithm = in_algorithm(0);
    accept_fraction = in_frac(0);
    batch_size = in_batch(0);
    epochs = in_epochs(0);
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
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>();
}


