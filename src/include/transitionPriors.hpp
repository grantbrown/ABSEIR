#ifndef SPATIALSEIR_TRANSITION_PRIORS
#define SPATIALSEIR_TRANSITION_PRIORS
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(transitionPriors)
class transitionPriors : public modelComponent
{
    public:
        transitionPriors(SEXP mode);
        void setPriorsFromProbabilities(SEXP p_ei, SEXP p_ir, SEXP p_ei_ess, SEXP p_ir_ess);
        void setPathSpecificPriors(SEXP Zmat1, SEXP Zmat2);
        void setUniformExpPriors();
        int getModelComponentType();
        void setPriorsManually(SEXP priorAlpha_gammaEI, SEXP priorBeta_gammaEI,
                               SEXP priorAlpha_gammaIR, SEXP priorBeta_gammaIR);
        void summary();

        Eigen::VectorXd gamma_ei_params;
        Eigen::VectorXd gamma_ir_params;

        Eigen::MatrixXd gamma_ei_ps_probs;
        Eigen::MatrixXd gamma_ir_ps_probs;

        ~transitionPriors();
};


#endif
