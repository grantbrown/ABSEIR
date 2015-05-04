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
        std::vector<int> S0;
        std::vector<int> E0;
        std::vector<int> I0;
        std::vector<int> R0;
        ~initialValueContainer();
};


#endif
