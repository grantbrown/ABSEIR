#ifndef SPATIALSEIR_SAMPLING_CONTROL
#define SPATIALSEIR_SAMPLING_CONTROL
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(samplingControl)
class samplingControl : public modelComponent
{
    public:
        samplingControl(SEXP in_width, SEXP in_seed, SEXP in_cores);
        ~samplingControl();
    int getModelComponentType();
    int simulation_width;
    int random_seed;
    int CPU_cores;
};


#endif
