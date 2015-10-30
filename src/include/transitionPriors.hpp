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
        void setUniformExpPriors();
        int getModelComponentType();
        void summary();

        Eigen::MatrixXd gamma_ei_params;
        Eigen::MatrixXd gamma_ir_params;

        std::string mode;
        int max_latent;
        int max_infectious;
        double inf_mean;

        ~transitionPriors();
};


#endif
