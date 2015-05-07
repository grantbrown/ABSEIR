#ifndef SPATIALSEIR_DISTANCE_MODEL
#define SPATIALSEIR_DISTANCE_MODEL
#include <Rcpp.h>
#include<modelComponent.hpp>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(distanceModel)
class distanceModel : public modelComponent
{
    public:
        distanceModel();
        virtual void addDistanceMatrix(NumericMatrix distMat);
        int getModelComponentType();
        virtual void summary();
        virtual void setPriorParameters(double alpha, double beta);
        virtual int getNumDistanceMatrices();

        int numLocations;
        Eigen::VectorXd spatial_prior;
        std::vector<Eigen::MatrixXd> dm_list;

        ~distanceModel();
};


#endif
