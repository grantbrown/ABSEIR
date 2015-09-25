#include <Rcpp.h>
#include <samplingControl.hpp>


using namespace Rcpp;

samplingControl::samplingControl(SEXP integerParameters,
                                 SEXP numericParameters)
{
    Rcpp::IntegerVector inIntegerParams(integerParameters);
    Rcpp::NumericVector inNumericParams(numericParameters);

    if (inIntegerParams.size() != 8 ||
        inNumericParams.size() != 2)
    {
        Rcpp::stop("Exactly 10 samplingControl parameters are required.");
    }

    simulation_width = inIntegerParams(0);
    random_seed = inIntegerParams(1); 
    CPU_cores = inIntegerParams(2); 
    algorithm = inIntegerParams(3);
    batch_size = inIntegerParams(4);
    epochs = inIntegerParams(5);
    max_batches = inIntegerParams(6);
    multivariatePerturbation = inIntegerParams(7) != 0;

    accept_fraction = inNumericParams(0);
    shrinkage = inNumericParams(1);

    if (algorithm != ALG_BasicABC && algorithm != ALG_ModifiedBeaumont2009)
    {
        Rcpp::stop("Algorithm specification must be of length 1 and equal to 1 or 2.");
    }
    if (max_batches <= 0)
    {
        Rcpp::stop("max_batches must be greater than zero.");
    }
}

samplingControl::~samplingControl()
{
    if (prot !=0 ){
        Rcpp::stop("can't delete samplingControl, still being used.\n");
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
    .constructor<SEXP, SEXP>();
}


