#ifndef SPATIALSEIR_MODEL_INC
#define SPATIALSEIR_MODEL_INC
#include <Rcpp.h>
#include "caf/all.hpp"
#include "./dataModel.hpp"
#include "./distanceModel.hpp"
#include "./exposureModel.hpp"
#include "./initialValueContainer.hpp"
#include "./reinfectionModel.hpp"
#include "./samplingControl.hpp"
#include "./spatialSEIRModel.hpp"
#include "./transitionPriors.hpp"

struct simulationResultSet
{
    Eigen::MatrixXi S;
    Eigen::MatrixXi E;
    Eigen::MatrixXi I;
    Eigen::MatrixXi R;
    Eigen::MatrixXi S_star;
    Eigen::MatrixXi E_star;
    Eigen::MatrixXi I_star;
    Eigen::MatrixXi R_star;
    Eigen::MatrixXd X;
    Eigen::MatrixXd p_se;
    Eigen::MatrixXd p_ei;
    Eigen::MatrixXd p_ir;
    Eigen::MatrixXd rho;
    Eigen::MatrixXd beta;
    double result;
};

class dataModel;
class exposureModel;
class distanceModel;
class initialValueContainer;
class reinfectionModel;
class samplingControl;
class transitionPriors;

using namespace Rcpp;
using namespace caf;

class spatialSEIRModel
{
    public: 
        //Constructor
        spatialSEIRModel(dataModel& dataModel_,
                         exposureModel& exposureModel_,
                         reinfectionModel& reinfectionModel_,
                         distanceModel& distanceModel_,
                         transitionPriors& transitionPriors_,
                         initialValueContainer& initialValueContainer_,
                         samplingControl& samplingControl_);

        Rcpp::NumericVector marginalPosteriorEstimates(SEXP inParams);
        Rcpp::List simulate(SEXP inParams);
        //Destructor
        ~spatialSEIRModel();
    private:
        int ncalls;
        dataModel* dataModelInstance;
        exposureModel* exposureModelInstance;
        reinfectionModel* reinfectionModelInstance;
        distanceModel* distanceModelInstance;
        transitionPriors* transitionPriorsInstance;
        initialValueContainer* initialValueContainerInstance;
        samplingControl* samplingControlInstance;
        scoped_actor* self;
};

#endif
