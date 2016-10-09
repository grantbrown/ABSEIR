#ifndef SPATIALSEIR_SAMPLING_CONTROL
#define SPATIALSEIR_SAMPLING_CONTROL

//#define SPATIALSEIR_SINGLETHREAD

#define ALG_BasicABC 1
#define ALG_ModifiedBeaumont2009 2
#define ALG_DelMoral2012 3
#define ALG_Simulate 4

#include <Rcpp.h>
#include<modelComponent.hpp>


using namespace Rcpp;

RCPP_EXPOSED_CLASS(samplingControl)
class samplingControl : public modelComponent
{
    public:
        samplingControl(SEXP integerParameters,
                        SEXP numericParameters);
        ~samplingControl();
    void summary();
    int getModelComponentType();
    int simulation_width;
    int random_seed;
    int algorithm;
    double target_eps;
    double accept_fraction;
    double shrinkage;
    int batch_size;
    int max_batches;
    int CPU_cores;
    int epochs;
    int m;
    bool multivariatePerturbation;
};


#endif
