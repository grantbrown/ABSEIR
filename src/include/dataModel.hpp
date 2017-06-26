#ifndef SPATIALSEIR_DATA_MODEL
#define SPATIALSEIR_DATA_MODEL
#include <Rcpp.h>
#include <modelComponent.hpp>
#include <Eigen/Core>

#define COMPARTMENT_I_STAR 0
#define COMPARTMENT_R_STAR 1
#define COMPARTMENT_I 2



using namespace Rcpp;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
RCPP_EXPOSED_CLASS(dataModel)
class dataModel : public modelComponent
{
    public:
        dataModel(SEXP Y, SEXP type, SEXP compartment, SEXP cumulative, 
                  SEXP phi, SEXP na_mask);
        virtual void summary();
        Rcpp::IntegerVector* compartmentDimensions;
        int getModelComponentType();
        int nLoc;
        int nTpt;
        int dataModelType;
        int dataModelCompartment;
        bool cumulative;
        double phi;
        double report_fraction;
        double report_fraction_ess;
        Eigen::MatrixXi Y;
        MatrixXb na_mask;
        ~dataModel();
};


#endif
