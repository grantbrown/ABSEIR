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
        void setPathSpecificPriors(SEXP Zmat1, SEXP Zmat2, SEXP inf_mean);
        void setPriorsForWeibull(SEXP E_to_I_params,
                                 SEXP I_to_R_params,
                                 SEXP maxEI,
                                 SEXP maxIR);
        void setUniformExpPriors();
        int getModelComponentType();
        void summary();

        Eigen::MatrixXd E_to_I_params;
        Eigen::MatrixXd I_to_R_params;

        std::string mode;
        int max_latent;
        int max_infectious;
        double inf_mean;

        ~transitionPriors();
};


#endif
