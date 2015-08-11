#ifndef SPATIALSEIR_SAMPLING_CONTROL
#define SPATIALSEIR_SAMPLING_CONTROL

#define ALG_BasicABC 1
#define ALG_ModifiedBeaumont2009 2

#include <Rcpp.h>
#include<modelComponent.hpp>


using namespace Rcpp;

RCPP_EXPOSED_CLASS(samplingControl)
class samplingControl : public modelComponent
{
    public:
        samplingControl(SEXP in_width, SEXP in_seed, SEXP in_cores,
                        SEXP in_fraction, SEXP in_batch_size,
                        SEXP in_algorithm, SEXP in_epochs);
        ~samplingControl();
    int getModelComponentType();
    int simulation_width;
    int random_seed;
    int algorithm;
    double accept_fraction;
    int batch_size;
    int CPU_cores;
    int epochs;
};


#endif
