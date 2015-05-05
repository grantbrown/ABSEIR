#ifndef SPATIALSEIR_TRANSITION_PRIORS
#define SPATIALSEIR_TRANSITION_PRIORS
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(transitionPriors)
class transitionPriors : public modelComponent
{
    public:
        transitionPriors();
        void setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, SEXP p_ei_ess, SEXP p_ir_ess);
        void setUniformPriors();
        int getModelComponentType();
        void setPriorsManually(SEXP priorAlpha_gammaEI, SEXP priorBeta_gammaEI,
                               SEXP priorAlpha_gammaIR, SEXP priorBeta_gammaIR);
        void summary();

        std::vector<double> gamma_ei_params;
        std::vector<double> gamma_ir_params;
        ~transitionPriors();
};


#endif
