#include <Rcpp.h>
#include <samplingControl.hpp>


using namespace Rcpp;

samplingControl::samplingControl(SEXP integerParameters,
                                 SEXP numericParameters)
{
    Rcpp::IntegerVector inIntegerParams(integerParameters);
    Rcpp::NumericVector inNumericParams(numericParameters);

    if (inIntegerParams.size() != 10 ||
        inNumericParams.size() != 4)
    {
        Rcpp::stop("Exactly 15 samplingControl parameters are required.");
    }

    simulation_width = inIntegerParams(0);
    random_seed = inIntegerParams(1); 
    CPU_cores = inIntegerParams(2); 
    algorithm = inIntegerParams(3);
    batch_size = inIntegerParams(4);
    init_batch_size = inIntegerParams(5);
    epochs = inIntegerParams(6);
    max_batches = inIntegerParams(7);
    multivariatePerturbation = inIntegerParams(8) != 0;
    m = inIntegerParams(9);
#ifdef SPATIALSEIR_SINGLETHREAD
    if (CPU_cores > 1)
    {
        Rcpp::stop("Error: multiple cores requested for ABSEIR compiled in single thread mode");
    }
#endif


    accept_fraction = inNumericParams(0);
    shrinkage = inNumericParams(1);
    lpow = inNumericParams(2);
    target_eps = inNumericParams(3);
    

    if (algorithm != ALG_BasicABC && 
        algorithm != ALG_ModifiedBeaumont2009 && 
        algorithm != ALG_DelMoral2012 && 
        algorithm != ALG_Simulate)
    {
        Rcpp::stop("Algorithm specification must be of length 1 and equal to 1 or 2 or 3.");
    }
    if (max_batches <= 0)
    {
        Rcpp::stop("max_batches must be greater than zero.");
    }
}

samplingControl::samplingControl(samplingControl* tocopy, int seed_offs)
{
    simulation_width = tocopy -> simulation_width;
    random_seed = (tocopy -> random_seed) + seed_offs;
    algorithm = (tocopy -> algorithm);
    target_eps = (tocopy -> target_eps);
    accept_fraction = (tocopy -> accept_fraction);
    shrinkage = (tocopy -> shrinkage);
    batch_size = (tocopy -> batch_size);
    init_batch_size = (tocopy -> init_batch_size); 
    max_batches = (tocopy -> max_batches);
    CPU_cores = (tocopy -> CPU_cores);
    epochs = (tocopy -> epochs);
    m = (tocopy -> m);
    lpow = (tocopy -> lpow);
    multivariatePerturbation = (tocopy -> multivariatePerturbation);
}

void samplingControl::summary()
{
    Rcpp::Rcout << "Sampling Control Summary\n";
    Rcpp::Rcout << "------------------------\n";
    Rcpp::Rcout << "    algorithm: " << algorithm << "\n";
    Rcpp::Rcout << "    simulation_width: " << simulation_width << "\n";
    Rcpp::Rcout << "    random_seed: " << random_seed << "\n";
    Rcpp::Rcout << "    CPU_cores: " << CPU_cores << "\n";
    Rcpp::Rcout << "    init_batch_size: " << init_batch_size << "\n";
    Rcpp::Rcout << "    batch_size: " << batch_size << "\n";
    Rcpp::Rcout << "    epochs: " << epochs << "\n";
    Rcpp::Rcout << "    max_batches: " << max_batches << "\n";
    Rcpp::Rcout << "    multivariatePerturbation: " << multivariatePerturbation << "\n";
    Rcpp::Rcout << "    m: " << m << "\n";
    Rcpp::Rcout << "    accept_fraction: " << accept_fraction << "\n";
    Rcpp::Rcout << "    shrinkage: " << shrinkage << "\n";
    Rcpp::Rcout << "    lpow: " << lpow << "\n";
    Rcpp::Rcout << "    target_eps: " << target_eps << "\n";
    Rcpp::Rcout << "    Note: not all parameters are used for all algorithms.\n\n";

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


