#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include <Eigen/Core>
#include <modelComponent.hpp>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(dataModel)
class dataModel : public modelComponent
{
    public:
        dataModel(SEXP Y, SEXP type, SEXP compartment);
        virtual void summary();
        virtual void setOverdispersionParameters(SEXP priorAlpha, SEXP priorBeta);
        Rcpp::IntegerVector* compartmentDimensions;
        int getModelComponentType();
        std::vector<double> priorParameters;
        std::vector<double> initialParameterValues;
        int nLoc;
        int nTpt;
        int dataModelType;
        int dataModelCompartment;
        Eigen::MatrixXi Y;
        ~dataModel();
};


#endif
