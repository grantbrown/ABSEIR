#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include <modelComponent.hpp>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(dataModel)
class dataModel : public modelComponent
{
    public:
        dataModel(SEXP Y, SEXP type, SEXP compartment, SEXP phi);
        virtual void summary();
        Rcpp::IntegerVector* compartmentDimensions;
        int getModelComponentType();
        int nLoc;
        int nTpt;
        int dataModelType;
        int dataModelCompartment;
        double phi;
        Eigen::MatrixXi Y;
        ~dataModel();
};


#endif
