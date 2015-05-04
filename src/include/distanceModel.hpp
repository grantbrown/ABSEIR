#ifndef SPATIALSEIR_DISTANCE_MODEL
#define SPATIALSEIR_DISTANCE_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>
#include <Eigen/Core>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(distanceModel)
class distanceModel : public modelComponent
{
    public:
        distanceModel();
        virtual void addDistanceMatrix(NumericMatrix distMat);
        int getModelComponentType();
        virtual void summary();
        virtual int getNumDistanceMatrices();

        int numLocations;
        double priorAlpha;
        double priorBeta;
        std::vector<Eigen::MatrixXd> dm_list;

        ~distanceModel();
};


#endif
