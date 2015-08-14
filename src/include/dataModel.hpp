#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include <modelComponent.hpp>
#include <Eigen/Core>

using namespace Rcpp;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
RCPP_EXPOSED_CLASS(dataModel)
class dataModel : public modelComponent
{
    public:
        dataModel(SEXP Y, SEXP type, SEXP compartment, SEXP phi,
                  SEXP na_mask);
        virtual void summary();
        Rcpp::IntegerVector* compartmentDimensions;
        int getModelComponentType();
        int nLoc;
        int nTpt;
        int dataModelType;
        int dataModelCompartment;
        double phi;
        Eigen::MatrixXi Y;
        MatrixXb na_mask;
        ~dataModel();
};


#endif
