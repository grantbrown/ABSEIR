#ifndef SPATIALSEIR_EXPOSURE_MODEL
#define SPATIALSEIR_EXPOSURE_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>
#include<Eigen/Core>

using namespace Rcpp;


RCPP_EXPOSED_CLASS(exposureModel)
class exposureModel : public modelComponent
{
    public:
        exposureModel(SEXP X, SEXP ntpt, SEXP nloc,SEXP priorMean, SEXP precision);
        int getModelComponentType();
        virtual Rcpp::NumericVector getOffset();
        virtual void setOffset(Rcpp::NumericVector offs);
        std::vector<double> offset;
        int nTpt;
        int nLoc;
        Eigen::MatrixXd X;
        std::vector<double> betaPriorPrecision;
        std::vector<double> betaPriorMean;
        ~exposureModel();
};


#endif
