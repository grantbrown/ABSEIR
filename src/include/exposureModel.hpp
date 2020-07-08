#ifndef SPATIALSEIR_EXPOSURE_MODEL
#define SPATIALSEIR_EXPOSURE_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;


RCPP_EXPOSED_CLASS(exposureModel)
class exposureModel : public modelComponent
{
    public:
        exposureModel(SEXP X, SEXP ntpt, SEXP nloc,SEXP priorMean, SEXP precision);
        exposureModel(exposureModel* tocopy);
        int getModelComponentType();
        virtual Rcpp::NumericVector getOffset();
        virtual void setOffset(Rcpp::NumericVector offs);
        virtual void summary();
        int nTpt;
        int nLoc;
        Eigen::MatrixXd X;
        Eigen::VectorXd betaPriorPrecision;
        Eigen::VectorXd betaPriorMean;
        Eigen::VectorXd offset;
        ~exposureModel();
};


#endif
