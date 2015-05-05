#ifndef SPATIALSEIR_INITVALUE_CONTAINER
#define SPATIALSEIR_INITVALUE_CONTAINER
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(initialValueContainer)

class initialValueContainer : public modelComponent
{
    public:
        initialValueContainer();
        void setInitialValues(SEXP S0, SEXP E0, SEXP I0, SEXP R0);
        int getModelComponentType();
        Eigen::VectorXi S0;
        Eigen::VectorXi E0;
        Eigen::VectorXi I0;
        Eigen::VectorXi R0;
        ~initialValueContainer();
};


#endif
