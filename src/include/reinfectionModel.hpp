#ifndef SPATIALSEIR_REINFECTION_MODEL
#define SPATIALSEIR_REINFECTION_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;


RCPP_EXPOSED_CLASS(reinfectionModel)
class reinfectionModel : public modelComponent
{
    public:
        reinfectionModel(SEXP reinfectionMode);
        reinfectionModel(reinfectionModel* tocopy);
        int getModelComponentType();
        virtual void summary();
        virtual void buildReinfectionModel(SEXP _X, SEXP _priorMean, SEXP _prec);
        int reinfectionMode;
        Eigen::MatrixXd X_rs;
        Eigen::VectorXd betaPriorPrecision;
        Eigen::VectorXd betaPriorMean;
        ~reinfectionModel();
};


#endif
