#ifndef SPATIALSEIR_REINFECTION_MODEL
#define SPATIALSEIR_REINFECTION_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>
#include<Eigen/Core>

using namespace Rcpp;


RCPP_EXPOSED_CLASS(reinfectionModel)
class reinfectionModel : public modelComponent
{
    public:
        reinfectionModel(SEXP reinfectionMode);
        int getModelComponentType();
        virtual void buildReinfectionModel(SEXP _X, SEXP _priorMean, SEXP _prec);
        int reinfectionMode;
        Eigen::MatrixXd X_rs;
        std::vector<double> betaPriorPrecision;
        std::vector<double> betaPriorMean;
        ~reinfectionModel();
};


#endif
