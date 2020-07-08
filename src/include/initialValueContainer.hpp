#ifndef SPATIALSEIR_INITVALUE_CONTAINER
#define SPATIALSEIR_INITVALUE_CONTAINER
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(initialValueContainer)

class initialValueContainer : public modelComponent
{
    public:
        initialValueContainer(int ivc_type);
        initialValueContainer(initialValueContainer* tocopy);

        void setInitialValues(SEXP S0, SEXP E0, SEXP I0, SEXP R0,
                              SEXP MS0_, SEXP ME0_, SEXP MIO_, SEXP MR0_);

        virtual void summary();
        int getModelComponentType();
        int type;
        Eigen::VectorXi S0;
        Eigen::VectorXi E0;
        Eigen::VectorXi I0;
        Eigen::VectorXi R0;
		Eigen::VectorXi S0_max;
        Eigen::VectorXi E0_max;
        Eigen::VectorXi I0_max;
        Eigen::VectorXi R0_max;
        
        ~initialValueContainer();
};


#endif
